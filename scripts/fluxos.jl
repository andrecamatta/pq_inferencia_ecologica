# Fluxo de votos 1º → 2º turno por seção eleitoral, com faixas de incerteza e saldos.
#
# Entrada: data/secoes/<ano>.arrow (scripts/preparar.jl).
# Saída: resultados/fluxos-secao-<ano>.json (mesmo formato lido pelo painel da apuração).
#
# Todas as versões usam efeito fixo de local de votação (comparam seções do mesmo local; a versão
# "fe_zona" compara seções da mesma zona). Faixa de cada fluxo: envelope das estimativas das
# versões, alargado pela distância entre a estimativa principal e os percentis 5 e 95 do seu
# bootstrap (municípios sorteados com reposição dentro de cada UF, levados às taxas),
#
#     bruta = [min_s N_s − (N_p − Q_p(5%)),  max_s N_s + (Q_p(95%) − N_p)],
#
# e calibrada pela validação por simulação (scripts/validar.jl, campo "calibracao"):
#
#     faixa = [N_p − c·(N_p − mín) − d·N_p,  N_p + c·(máx − N_p) + d·N_p].
#
# Saldos: para cada par de grupos presentes nos dois turnos, S = N(A→B) − N(B→A), com faixa
# montada do mesmo jeito (bootstrap do próprio saldo, d aplicado à soma ida + volta).
#
# Uso: julia -t 8 --project=. scripts/fluxos.jl 2026 [reamostragens]

using Arrow, DataFrames, JSON3, Statistics, Dates, Printf, LinearAlgebra
# O paralelismo é feito pelas threads do Julia; o BLAS em várias threads ao mesmo tempo só disputa CPU.
Threads.nthreads() > 1 && BLAS.set_num_threads(1)
include(joinpath(dirname(@__DIR__), "src", "InferenciaEcologica.jl"))
using .InferenciaEcologica
const IE = InferenciaEcologica

const RAIZ = dirname(@__DIR__)
const XC = [:x_lula, :x_outros, :x_adv, :x_bn, :x_abst]
const YC = [:y_lula, :y_adv, :y_bn, :y_abst]
const SOCIAIS = [:esc_baixa, :esc_superior, :idade_60]
const ADVERSARIO = Dict("2022" => "Bolsonaro", "2026" => "Flávio")
const PRINCIPAL = "fe_social"
const ESPECS = [  # nome => (agrupamento, covariáveis, efeito fixo, descrição)
    "fe_social" => (:uf, SOCIAIS, :local, "Por estado + escolaridade e idade, comparando seções do mesmo local (principal)"),
    "fe_uf" => (:uf, Symbol[], :local, "Por estado, comparando seções do mesmo local"),
    "fe_regiao" => (:regiao, Symbol[], :local, "Por região, comparando seções do mesmo local"),
    "fe_nacional" => (nothing, Symbol[], :local, "Uma matriz para o país, comparando seções do mesmo local"),
    "fe_zona" => (:uf, SOCIAIS, :zona, "Por estado + escolaridade e idade, comparando seções da mesma zona"),
]
const DESTINO = Dict(1 => 1, 3 => 2, 4 => 3, 5 => 4)           # origem que também existe no 2º turno → coluna
const PARES_SALDO = [(3, 1), (5, 1), (5, 3), (4, 1), (4, 3), (5, 4)]   # saldo a favor do segundo (índices 1-based)

fatoriza(v) = (m = Dict{eltype(v),Int}(); [get!(m, x, length(m) + 1) for x in v])

function carregar(ano)
    s = DataFrame(Arrow.Table(get(ENV, "SECOES", joinpath(RAIZ, "data", "secoes", "$(ano).arrow"))))
    X, Y = Float64.(Matrix(s[!, XC])), Float64.(Matrix(s[!, YC]))
    n = vec(sum(X; dims=2))
    cov = Dict{Symbol,Vector{Float64}}(c => Float64.(s[!, c]) for c in SOCIAIS)
    cov[:forca] = X[:, 1] ./ max.(X[:, 1] .+ X[:, 3], 1)
    cov[:bn] = X[:, 4] ./ n
    for (c, v) in cov
        cov[c] = v .- sum(v .* n) / sum(n)
    end
    uf = String.(s.uf)
    cluster = fatoriza(uf .* string.(s.municipio))
    # Locais em que alguma seção mudou de eleitorado entre os turnos ficam fora das diferenças.
    chave = uf .* string.(s.municipio) .* "_" .* string.(s.zona)
    localid = fatoriza(chave .* "_" .* string.(s.local_votacao))
    mudou = falses(maximum(localid))
    for i in eachindex(localid)
        s.aptos1[i] != s.aptos2[i] && (mudou[localid[i]] = true)
    end
    local_ = [mudou[localid[i]] || s.local_votacao[i] < 0 ? -1 : localid[i] for i in eachindex(localid)]
    (; X, Y, uf, cov, cluster, efeito=Dict(:local => local_, :zona => fatoriza(chave)), localid)
end

function ajustar(D, espec; idx=nothing, Y=D.Y)
    grupo, covs, fe, _ = Dict(ESPECS)[espec]
    X, uf = D.X, D.uf
    Z = isempty(covs) ? nothing : hcat((D.cov[c] for c in covs)...)
    loc = D.efeito[fe]
    if idx !== nothing
        X, Y, uf, loc = X[idx, :], Y[idx, :], uf[idx], loc[idx]
        Z = Z === nothing ? nothing : Z[idx, :]
    end
    agrega(fluxos(X, Y, uf, Z, grupo; local_=loc), uf)
end

saldo(M::AbstractMatrix, a, b) = M[a, DESTINO[b]] - M[b, DESTINO[a]]

"Calibração (c, d) gravada pela validação; sem ela, a faixa sai bruta."
function calibracao()
    arq = get(ENV, "CALIBRACAO", joinpath(RAIZ, "resultados", "validacao.json"))
    isfile(arq) || return nothing
    cal = JSON3.read(read(arq, String))
    haskey(cal, :calibracao) ? cal.calibracao : nothing
end

function main(ano, nboot)
    t0 = time()
    D = carregar(ano)
    println("$(size(D.X, 1)) seções, $(Threads.nthreads()) threads"); flush(stdout)
    nomes = first.(ESPECS)
    res = Vector{Dict{String,Matrix{Float64}}}(undef, length(nomes))
    Threads.@threads for e in eachindex(nomes)
        res[e] = ajustar(D, nomes[e])
    end
    pontos = Dict(zip(nomes, res))
    @printf("especificações ajustadas (%.0fs)\n", time() - t0)
    boots = Vector{Dict{String,Matrix{Float64}}}(undef, nboot)
    feitos = Threads.Atomic{Int}(0)
    Threads.@threads for b in 1:nboot
        rs = Reamostrador(D.uf, D.cluster, b)
        boots[b] = na_escala(ajustar(D, PRINCIPAL; idx=amostra(rs)), pontos[PRINCIPAL])
        Threads.atomic_add!(feitos, 1) % 50 == 49 && @printf("  %d/%d reamostragens (%.0fs)\n", feitos[], nboot, time() - t0)
    end
    @printf("bootstrap: %d reamostragens (%.0fs)\n", nboot, time() - t0)

    cal = calibracao()
    c, d = cal === nothing ? (1.0, 0.0) : (Float64(cal.c), Float64(cal.d))
    r(M) = [[round(Int, M[i, j]) for j in axes(M, 2)] for i in axes(M, 1)]
    niveis = Dict{String,Any}()
    for (nome, N) in pontos[PRINCIPAL]
        B = [b[nome] for b in boots]
        saldos = Dict{String,Any}()
        for (a, b) in PARES_SALDO
            S, Sb = saldo(N, a, b), [saldo(M, a, b) for M in B]
            s05, s95 = quantile(Sb, [0.05, 0.95])
            ss = [saldo(pontos[e][nome], a, b) for e in nomes]
            lo, hi = min(minimum(ss) - max(S - s05, 0), S), max(maximum(ss) + max(s95 - S, 0), S)
            vol = N[a, DESTINO[b]] + N[b, DESTINO[a]]
            saldos["$(a-1)-$(b-1)"] = Dict("saldo" => round(Int, S), "faixa" => [round(Int, S - c * (S - lo) - d * vol), round(Int, S + c * (hi - S) + d * vol)],
                                           "ida" => round(Int, N[a, DESTINO[b]]), "volta" => round(Int, N[b, DESTINO[a]]))
        end
        q05 = [quantile([M[i, j] for M in B], 0.05) for i in 1:5, j in 1:4]
        q95 = [quantile([M[i, j] for M in B], 0.95) for i in 1:5, j in 1:4]
        mn = reduce((A, e) -> min.(A, pontos[e][nome]), nomes; init=fill(Inf, 5, 4))
        mx = reduce((A, e) -> max.(A, pontos[e][nome]), nomes; init=fill(-Inf, 5, 4))
        lo = min.(max.(mn .- max.(N .- q05, 0), 0), N)
        hi = max.(mx .+ max.(q95 .- N, 0), N)
        niveis[nome] = Dict(
            "fluxos" => r(N), "faixa" => [r(max.(N .- c .* (N .- lo) .- d .* N, 0)), r(N .+ c .* (hi .- N) .+ d .* N)],
            "faixaBruta" => [r(lo), r(hi)], "amostragem" => [r(q05), r(q95)],
            "porEspec" => Dict(e => r(pontos[e][nome]) for e in nomes),
            "origem" => round.(Int, vec(sum(N; dims=2))), "destino" => round.(Int, vec(sum(N; dims=1))), "saldos" => saldos)
    end
    adv = get(ADVERSARIO, ano, "Adversário")
    out = Dict("geradoEm" => string(now(UTC)) * "Z", "ano" => parse(Int, ano), "unidade" => "seção",
               "secoes" => size(D.X, 1), "bootstrap" => nboot, "bootstrapAlt" => 0, "kUf" => IE.K_GRUPO,
               "origens" => ["Lula", "Outros candidatos", adv, "Brancos e nulos", "Abstenção"],
               "destinos" => ["Lula", adv, "Brancos e nulos", "Abstenção"], "calibracao" => cal,
               "especs" => Dict(e => last(v) for (e, v) in ESPECS), "principal" => PRINCIPAL, "niveis" => niveis,
               "implementacao" => "julia")
    mkpath(joinpath(RAIZ, "resultados"))
    write(joinpath(RAIZ, "resultados", "fluxos-secao-$(ano).json"), JSON3.write(out))
    br = niveis["Brasil"]
    @printf("Brasil (milhões) [faixa] (%.0fs):\n", time() - t0)
    for (i, o) in enumerate(out["origens"])
        println("  ", rpad(o, 18), join([@sprintf("%-6s %5.2f [%4.2f–%4.2f]", first(out["destinos"][j], 6), br["fluxos"][i][j] / 1e6,
                                                  br["faixa"][1][i][j] / 1e6, br["faixa"][2][i][j] / 1e6) for j in 1:4], "  "))
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__()) && main(ARGS[1], length(ARGS) > 1 ? parse(Int, ARGS[2]) : 400)
