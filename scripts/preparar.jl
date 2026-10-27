# Monta a tabela por seção (1º e 2º turno de presidente + perfil do eleitorado) de um ano.
#
# Entrada: arquivos do TSE em data/brutos/<ano>/ (scripts/baixar.jl), ou na pasta indicada pela
# variável de ambiente TSE_BRUTOS.
# Saída: data/secoes/<ano>.arrow, uma linha por seção, com
#   x_lula, x_outros, x_adv, x_bn, x_abst   (1º turno)
#   y_lula, y_adv, y_bn, y_abst             (2º turno)
#   esc_baixa, esc_superior, idade_60       (frações do eleitorado da seção)
#   uf, municipio, zona, secao, local_votacao, area, aptos1, aptos2
#
# Uso: julia -t 8 --project=. scripts/preparar.jl 2026

using Arrow, DataFrames, ZipArchives, Mmap, Statistics

const RAIZ = dirname(@__DIR__)
const LULA, ADV = 13, 22
const Chave = Tuple{String,Int,Int,Int}   # uf, município, zona, seção

"Campos pedidos de uma linha do CSV: posição inicial e final de cada um no buffer, sem aspas."
struct Linha
    b::Vector{UInt8}
    ini::Vector{Int}
    fim::Vector{Int}
end

"Número inteiro do campo `k` (0 se vazio ou não numérico, como `pd.to_numeric(errors=\"coerce\")`)."
function num(l::Linha, k)
    v, neg = 0, false
    @inbounds for p in l.ini[k]:l.fim[k]
        c = l.b[p]
        if c == UInt8('-') && p == l.ini[k]
            neg = true
        elseif UInt8('0') <= c <= UInt8('9')
            v = 10v + (c - UInt8('0'))
        else
            return 0
        end
    end
    neg ? -v : v
end
txt(l::Linha, k) = String(l.b[l.ini[k]:l.fim[k]])

"""
Percorre um CSV `;`-separado dentro do zip, em blocos, chamando `f(linha)` para cada linha.
`nomes` são as colunas pedidas; na linha, a coluna `nomes[k]` é o campo `k`.
"""
function percorre(f, zipcaminho, filtro_nome, nomes)
    zr = ZipReader(mmap(zipcaminho))
    zip_openentry(zr, first(filter(filtro_nome, zip_names(zr)))) do arq
        cab = [strip(c, '"') for c in split(readline(arq), ';')]
        pos = Dict(String(c) => i for (i, c) in enumerate(cab))
        alvo = zeros(Int, length(cab))
        for (k, n) in enumerate(nomes)
            alvo[pos[n]] = k
        end
        l = Linha(UInt8[], zeros(Int, length(nomes)), zeros(Int, length(nomes)))
        resto = UInt8[]
        bloco = Vector{UInt8}(undef, 1 << 25)
        while !eof(arq)
            n = readbytes!(arq, bloco)
            dados = vcat(resto, view(bloco, 1:n))
            ultimo = something(findlast(==(UInt8('\n')), dados), 0)
            if eof(arq) && ultimo < length(dados)   # última linha sem quebra
                push!(dados, UInt8('\n')); ultimo = length(dados)
            end
            l = Linha(dados, l.ini, l.fim)
            s = 1
            @inbounds while s <= ultimo
                if dados[s] == UInt8('\n') || dados[s] == UInt8('\r')   # linha vazia
                    s += 1; continue
                end
                e = s
                campo, ini = 1, s
                while true
                    c = dados[e]
                    if c == UInt8(';') || c == UInt8('\n') || c == UInt8('\r')
                        k = campo <= length(alvo) ? alvo[campo] : 0
                        if k > 0
                            a, z = ini, e - 1
                            a <= z && dados[a] == UInt8('"') && (a += 1)
                            a <= z && dados[z] == UInt8('"') && (z -= 1)
                            l.ini[k], l.fim[k] = a, z
                        end
                        c == UInt8(';') || break
                        campo += 1; ini = e + 1
                    end
                    e += 1
                end
                while e <= ultimo && dados[e] != UInt8('\n')
                    e += 1
                end
                f(l)
                s = e + 1
            end
            resto = dados[ultimo+1:end]
        end
    end
end

"Nomes das colunas do CSV dentro do zip."
function cabecalho(zipcaminho, filtro_nome)
    zr = ZipReader(mmap(zipcaminho))
    zip_openentry(io -> [String(strip(c, '"')) for c in split(readline(io), ';')], zr, first(filter(filtro_nome, zip_names(zr))))
end

const CHAVE_COLS = ["SG_UF", "CD_MUNICIPIO", "NR_ZONA", "NR_SECAO"]
chave(l, o=0) = (lowercase(txt(l, o + 1)), num(l, o + 2), num(l, o + 3), num(l, o + 4))

function votos(pasta, ano)
    v = Dict{Tuple{Int,Chave},NTuple{2,Int}}()   # (turno, chave) → (Lula, adversário)
    loc = Dict{Chave,Int}()                       # local de votação no 1º turno
    cols = vcat(CHAVE_COLS, ["CD_CARGO", "NR_TURNO", "NR_VOTAVEL", "QT_VOTOS", "NR_LOCAL_VOTACAO"])
    percorre(joinpath(pasta, "votacao_secao_$(ano)_BR.zip"), endswith(".csv"), cols) do l
        num(l, 5) == 1 || return
        turno, k, cand, q = num(l, 6), chave(l), num(l, 7), num(l, 8)
        a, b = get(v, (turno, k), (0, 0))
        v[(turno, k)] = (a + (cand == LULA) * q, b + (cand == ADV) * q)
        turno == 1 && !haskey(loc, k) && (loc[k] = num(l, 9))
    end
    v, loc
end

function detalhe(pasta, ano)
    d = Dict{Tuple{Int,Chave},NTuple{4,Int}}()   # aptos, abstenção, nominais, brancos + nulos + anulados
    cols = vcat(CHAVE_COLS, ["CD_CARGO", "NR_TURNO", "QT_APTOS", "QT_ABSTENCOES", "QT_VOTOS_NOMINAIS",
                             "QT_VOTOS_BRANCOS", "QT_VOTOS_NULOS", "QT_VOTOS_ANULADOS_APU_SEP"])
    # 2022: arquivo nacional *_BRASIL.csv; 2026: *_BR.csv (o zip de 2026 traz os dois; usa o primeiro, como o original)
    percorre(joinpath(pasta, "detalhe_votacao_secao_$(ano).zip"), n -> endswith(n, "_BRASIL.csv") || endswith(n, "_BR.csv"), cols) do l
        num(l, 5) == 1 || return
        d[(num(l, 6), chave(l))] = (num(l, 7), num(l, 8), num(l, 9), num(l, 10) + num(l, 11) + num(l, 12))
    end
    d
end

function perfil(pasta, ano)
    arquivos = sort(filter(f -> startswith(f, "perfil_eleitor_secao_$(ano)_") && endswith(f, ".zip"), readdir(pasta)))
    partes = Vector{Dict{Chave,NTuple{4,Int}}}(undef, length(arquivos))
    Threads.@threads for j in eachindex(arquivos)
        z = joinpath(pasta, arquivos[j])
        qt = "QT_ELEITORES_PERFIL" in cabecalho(z, endswith(".csv")) ? "QT_ELEITORES_PERFIL" : "QT_ELEITORES"   # 2022 / 2026
        p = Dict{Chave,NTuple{4,Int}}()   # n, baixa escolaridade, superior, 60+
        percorre(z, endswith(".csv"), vcat(CHAVE_COLS, [qt, "CD_GRAU_ESCOLARIDADE", "CD_FAIXA_ETARIA"])) do l
            k, n, esc, idade = chave(l), num(l, 5), num(l, 6), num(l, 7)
            a = get(p, k, (0, 0, 0, 0))
            p[k] = (a[1] + n, a[2] + n * (1 <= esc <= 3), a[3] + n * (7 <= esc <= 8), a[4] + n * (idade >= 6064))
        end
        partes[j] = p
        println("perfil ", arquivos[j][end-5:end-4])
    end
    merge(partes...)
end

media(v) = (w = filter(!isnan, v); isempty(w) ? NaN : mean(w))

function main(ano)
    pasta = get(ENV, "TSE_BRUTOS", joinpath(RAIZ, "data", "brutos", ano))
    (v, loc), d = votos(pasta, ano), detalhe(pasta, ano)
    println("votos: $(length(v)) seções×turno; detalhe: $(length(d))")
    linhas = NamedTuple[]
    for ((turno, k), r1) in d
        turno == 1 || continue
        r2 = get(d, (2, k), nothing)
        r2 === nothing && continue
        l1, a1 = get(v, (1, k), (0, 0))
        l2, a2 = get(v, (2, k), (0, 0))
        x = (l1, max(r1[3] - l1 - a1, 0), a1, r1[4], r1[2])
        y = (l2, a2, r2[4], r2[2])
        # Seções sem eleitores num dos turnos (agregadas a outra seção, por exemplo) ficam de fora.
        (sum(x) > 0 && sum(y) > 0) || continue
        push!(linhas, (uf=k[1], municipio=k[2], zona=k[3], secao=k[4], x_lula=x[1], x_outros=x[2], x_adv=x[3], x_bn=x[4],
                       x_abst=x[5], aptos1=r1[1], y_lula=y[1], y_adv=y[2], y_bn=y[3], y_abst=y[4], aptos2=r2[1],
                       local_votacao=get(loc, k, -1)))
    end
    s = sort!(DataFrame(linhas), [:uf, :municipio, :zona, :secao])
    pf = perfil(pasta, ano)
    prop(k, j) = (p = get(pf, k, nothing); p === nothing || p[1] == 0 ? NaN : p[j] / p[1])
    chaves = [(r.uf, r.municipio, r.zona, r.secao) for r in eachrow(s)]
    s.esc_baixa = prop.(chaves, 2); s.esc_superior = prop.(chaves, 3); s.idade_60 = prop.(chaves, 4)
    # Área como no painel: município, ou zona nas cidades com mais de uma zona.
    transform!(groupby(s, [:uf, :municipio]), :zona => (z -> length(unique(z))) => :nz)
    s.area = [r.uf * lpad(r.municipio, 5, '0') * (r.nz > 1 ? "z" * lpad(r.zona, 4, '0') : "") for r in eachrow(s)]
    select!(s, Not(:nz))
    sem_perfil = count(isnan, s.esc_baixa)
    for c in (:esc_baixa, :esc_superior, :idade_60)
        transform!(groupby(s, :area), c => (v -> (m = media(v); [isnan(x) ? m : x for x in v])) => c)
        geral = media(s[!, c])
        s[!, c] = [isnan(x) ? geral : x for x in s[!, c]]
    end
    mkpath(joinpath(RAIZ, "data", "secoes"))
    Arrow.write(joinpath(RAIZ, "data", "secoes", "$(ano).arrow"), s)
    X = [sum(s[!, c]) for c in (:x_lula, :x_outros, :x_adv, :x_bn, :x_abst)]
    Y = [sum(s[!, c]) for c in (:y_lula, :y_adv, :y_bn, :y_abst)]
    println("$(nrow(s)) seções ($sem_perfil sem perfil, preenchidas pela média da área); $(length(unique(s.area))) áreas")
    println("1º turno (mi): ", round.(X ./ 1e6; digits=2))
    println("2º turno (mi): ", round.(Y ./ 1e6; digits=2))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__()) && main(ARGS[1])
