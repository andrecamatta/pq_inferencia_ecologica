# Validação por simulação: com que frequência a faixa contém o fluxo verdadeiro?
#
# Usa a estrutura real das seções de 2022 (eleitorado do 1º turno, perfil, estados) e gera
# eleições sintéticas em que a transição de cada seção é conhecida. Em cada uma, roda o mesmo
# procedimento de scripts/fluxos.jl (versões + bootstrap por município) e compara a faixa com a
# verdade.
#
# Cenários:
#   A  modelo correto: a transição de cada seção é a estimada nos dados reais.
#   B  heterogeneidade não observada: desvios por município e por seção (escala logit).
#   C  efeito de contexto por seção: B + onde a SEÇÃO tem mais brancos e nulos, eleitores de
#      OUTROS candidatos e abstencionistas vão mais para Lula. Caso extremo: o contexto muda de
#      uma seção para outra da mesma escola.
#   D  efeito de contexto por local: o mesmo, ligado à taxa de brancos e nulos do LOCAL de votação.
#
# Calibração: a faixa publicada é alargada,
#     faixa = [N − c·(N − mín) − d·N,  N + c·(máx − N) + d·N],
# com (c, d) da menor largura que contém a verdade em 90% dos fluxos relevantes nos cenários A, B
# e D, ajustados nas réplicas pares e conferidos nas ímpares.
#
# Saída: resultados/validacao.json. Uso: julia -t 8 --project=. scripts/validar.jl [réplicas] [reamostragens]

using JSON3, Random, Statistics, Distributions, Dates, Printf
include(joinpath(@__DIR__, "fluxos.jl"))   # carrega o modelo, ESPECS, ajustar (sem executar main)

const FRACAO_MUN = 0.25
const CENARIOS = ["A" => "Modelo correto", "B" => "Heterogeneidade não observada (municípios e seções)",
                  "C" => "Efeito de contexto por seção (extremo)", "D" => "Efeito de contexto por local de votação"]
const SIG_MUN, SIG_SEC, CONTEXTO = 0.4, 0.3, 1.2
const RECORTES = ["Brasil", "Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul"]
const CALIBRA, ALVO = "ABD", 0.90

"Subamostra de municípios (por UF) e transições-base estimadas nos dados reais (sem efeito fixo)."
function base()
    d = carregar("2022")
    rng = MersenneTwister(7)
    manter = falses(size(d.X, 1))
    for u in unique(d.uf)
        s = d.uf .== u
        muns = unique(d.cluster[s])
        esc = Set(muns[randperm(rng, length(muns))[1:max(1, round(Int, length(muns) * FRACAO_MUN))]])
        manter .|= s .& in.(d.cluster, Ref(esc))
    end
    i = findall(manter)
    d = (; X=d.X[i, :], Y=d.Y[i, :], uf=d.uf[i], cov=Dict(c => v[i] for (c, v) in d.cov), cluster=fatoriza(d.cluster[i]),
         efeito=Dict(k => v[i] for (k, v) in d.efeito), localid=fatoriza(d.localid[i]))
    Z = hcat((d.cov[c] for c in SOCIAIS)...)
    _, M = fluxos(d.X, d.Y, d.uf, Z, :uf; devolver_matrizes=true)
    (; d..., M)
end

"Transições verdadeiras por seção e contagens do 2º turno sorteadas a partir delas."
function verdade(B, cenario, rng)
    M = copy(B.M)
    if cenario != "A"
        L = log.(M .+ 1e-4)
        mun = randn(rng, maximum(B.cluster), 5, 4) .* SIG_MUN
        L .+= mun[B.cluster, :, :] .+ randn(rng, size(L)) .* SIG_SEC
        if cenario in ("C", "D")
            zbn = copy(B.cov[:bn])
            if cenario == "D"   # taxa de brancos e nulos do local, não da seção
                n = vec(sum(B.X; dims=2))
                num, den = zeros(maximum(B.localid)), zeros(maximum(B.localid))
                for (k, l) in enumerate(B.localid)
                    num[l] += B.X[k, 4]; den[l] += n[k]
                end
                bnl = (num./den)[B.localid]
                zbn = bnl .- sum(bnl .* n) / sum(n)
            end
            zbn ./= std(zbn; corrected=false)
            for i in (2, 5)   # outros candidatos e abstenção → Lula
                L[:, i, 1] .+= CONTEXTO .* clamp.(zbn, -2, 3)
            end
        end
        M = exp.(L .- maximum(L; dims=3))
        M ./= sum(M; dims=3)
    end
    N = zeros(size(M))
    for a in axes(M, 1), i in 1:5
        x = round(Int, B.X[a, i])
        x > 0 && (N[a, i, :] = rand(rng, Multinomial(x, M[a, i, :] ./ sum(M[a, i, :]))))
    end
    dropdims(sum(N; dims=2); dims=2), N
end

function replica(B, cenario, k, nboot)
    rng = MersenneTwister(1000 * Int(cenario[1]) + k)
    Y, N = verdade(B, cenario, rng)
    real = agrega(N, B.uf)
    pontos = Dict(e => ajustar(B, e; Y) for e in first.(ESPECS))
    rs = Reamostrador(B.uf, B.cluster, k)
    boots = [na_escala(ajustar(B, PRINCIPAL; idx=amostra(rs), Y), pontos[PRINCIPAL]) for _ in 1:nboot]
    out = Dict{String,Any}()
    for nome in RECORTES
        P = pontos[PRINCIPAL][nome]
        q05 = [quantile([b[nome][i, j] for b in boots], 0.05) for i in 1:5, j in 1:4]
        q95 = [quantile([b[nome][i, j] for b in boots], 0.95) for i in 1:5, j in 1:4]
        mn = reduce((A, e) -> min.(A, pontos[e][nome]), first.(ESPECS); init=fill(Inf, 5, 4))
        mx = reduce((A, e) -> max.(A, pontos[e][nome]), first.(ESPECS); init=fill(-Inf, 5, 4))
        out[nome] = (verdade=real[nome], est=P, lo=min.(max.(mn .- max.(P .- q05, 0), 0), P), hi=max.(mx .+ max.(q95 .- P, 0), P),
                     q05, q95, porEspec=Dict(e => pontos[e][nome] for e in first.(ESPECS)))
    end
    (cenario, k, out)
end

alarga(E, LO, HI, c, d) = (max.(E .- c .* (E .- LO) .- d .* E, 0), E .+ c .* (HI .- E) .+ d .* E)

"Tensores réplicas × recortes × 5 × 4 de um campo."
tensor(rs, campo) = cat([cat([r[n][campo] for n in RECORTES]...; dims=3) for r in rs]...; dims=4)   # 5×4×recortes×réplicas

function calibrar(res)
    sel(par) = [r for (c, k, r) in res if occursin(c, CALIBRA) && k % 2 == par]
    V, E, LO, HI = (tensor(sel(0), f) for f in (:verdade, :est, :lo, :hi))
    tot = sum(V; dims=(1, 2))
    rel = V .>= 0.001 .* tot
    melhor = nothing
    for c in 1.0:0.1:5.0, d in 0.0:0.001:0.03
        lo, hi = alarga(E, LO, HI, c, d)
        mean(((V .>= lo) .& (V .<= hi))[rel]) >= ALVO || continue
        w = mean((hi .- lo) ./ tot)
        (melhor === nothing || w < melhor[1]) && (melhor = (w, round(c; digits=2), round(d; digits=4)))
    end
    c, d = melhor === nothing ? (5.0, 0.03) : melhor[2:3]
    V, E, LO, HI = (tensor(sel(1), f) for f in (:verdade, :est, :lo, :hi))
    lo, hi = alarga(E, LO, HI, c, d)
    rel = V .>= 0.001 .* sum(V; dims=(1, 2))
    Dict("c" => c, "d" => d, "alvo" => ALVO, "cenarios" => CALIBRA, "coberturaValidacao" => mean(((V .>= lo) .& (V .<= hi))[rel]))
end

function resumo(res, cal)
    sai = Dict{String,Any}()
    for (c, desc) in CENARIOS
        rs = [r for (cc, _, r) in res if cc == c]
        V, E, LO, HI, Q05, Q95 = (tensor(rs, f) for f in (:verdade, :est, :lo, :hi, :q05, :q95))
        tot = sum(V; dims=(1, 2))
        dentro, dentro_b = (V .>= LO) .& (V .<= HI), (V .>= Q05) .& (V .<= Q95)
        LOc, HIc = alarga(E, LO, HI, cal["c"], cal["d"])
        dentro_c = (V .>= LOc) .& (V .<= HIc)
        rel = V .>= 0.001 .* tot
        br(A) = A[:, :, 1, :]   # Brasil
        vies = dropdims(mean(br(E) .- br(V); dims=3); dims=3)
        pe(e) = cat([r["Brasil"].porEspec[e] for r in rs]...; dims=3)
        sai[c] = Dict("descricao" => desc, "replicas" => length(rs),
            "coberturaFaixa" => mean(dentro), "coberturaAmostragem" => mean(dentro_b),
            "coberturaFaixaRelevantes" => mean(dentro[rel]), "coberturaAmostragemRelevantes" => mean(dentro_b[rel]),
            "coberturaCalibradaRelevantes" => mean(dentro_c[rel]), "larguraCalibradaPctEleitorado" => mean(100 .* (HIc .- LOc) ./ tot),
            "porEspec" => Dict(e => Dict("erroMedioMil" => mean(abs.(pe(e) .- br(V))) / 1e3,
                                         "viesBnLulaMil" => mean(pe(e)[4, 1, :] .- br(V)[4, 1, :]) / 1e3) for e in first.(ESPECS)),
            "erroMedioPctEleitorado" => mean(100 .* abs.(E .- V) ./ tot), "larguraMediaPctEleitorado" => mean(100 .* (HI .- LO) ./ tot),
            "brasil" => Dict("cobertura" => [[round(mean(br(dentro)[i, j, :]); digits=3) for j in 1:4] for i in 1:5],
                             "coberturaCalibrada" => [[round(mean(br(dentro_c)[i, j, :]); digits=3) for j in 1:4] for i in 1:5],
                             "vies" => [[round(Int, vies[i, j]) for j in 1:4] for i in 1:5]))
    end
    sai
end

function principal(nrep, nboot)
    t0 = time()
    B = base()
    @printf("base: %d seções em %d municípios (%.0fs)\n", size(B.X, 1), maximum(B.cluster), time() - t0)
    tarefas = [(c, k) for (c, _) in CENARIOS for k in 0:nrep-1]
    res = Vector{Any}(undef, length(tarefas))
    feitos = Threads.Atomic{Int}(0)
    Threads.@threads for t in eachindex(tarefas)
        res[t] = replica(B, tarefas[t]..., nboot)
        Threads.atomic_add!(feitos, 1) % 8 == 7 && @printf("%d/%d réplicas (%.0fs)\n", feitos[], length(tarefas), time() - t0)
    end
    cal = calibrar(res)
    s = resumo(res, cal)
    out = Dict("geradoEm" => string(now(UTC)) * "Z", "replicas" => nrep, "bootstrap" => nboot, "secoes" => size(B.X, 1),
               "fracaoMunicipios" => FRACAO_MUN, "parametros" => Dict("sigmaMunicipio" => SIG_MUN, "sigmaSecao" => SIG_SEC, "contexto" => CONTEXTO),
               "calibracao" => cal, "cenarios" => s, "implementacao" => "julia")
    mkpath(joinpath(RAIZ, "resultados"))
    write(joinpath(RAIZ, "resultados", "validacao.json"), JSON3.write(out))
    for (c, _) in CENARIOS
        v = s[c]
        @printf("%s: cobertura relevantes %.0f%% (calibrada %.0f%%) | erro médio %.2f%% | viés BN→Lula %+.0f mil\n", c,
                100v["coberturaFaixaRelevantes"], 100v["coberturaCalibradaRelevantes"], v["erroMedioPctEleitorado"], v["brasil"]["vies"][4][1] / 1e3)
    end
    @printf("calibração: c = %.2f, d = %.4f; cobertura nas réplicas ímpares %.0f%%\ntotal %.0fs\n", cal["c"], cal["d"], 100cal["coberturaValidacao"], time() - t0)
end

principal(length(ARGS) > 0 ? parse(Int, ARGS[1]) : 24, length(ARGS) > 1 ? parse(Int, ARGS[2]) : 30)
