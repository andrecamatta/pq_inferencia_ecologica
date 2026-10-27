"""
Inferência ecológica do fluxo de votos entre turnos, por seção eleitoral.

Para a seção `a`, `x_a` são as contagens do 1º turno em 5 grupos (Lula, outros, adversário,
brancos e nulos, abstenção) e `y_a` as do 2º em 4 (Lula, adversário, brancos e nulos,
abstenção). O modelo é

    y_a / E_a ≈ B_aᵀ (x_a / E_a),   B_a = B⁽⁰⁾_g + Σ_k z_ak B⁽ᵏ⁾_g,

com linhas de B⁽⁰⁾ somando 1, linhas de B⁽ᵏ⁾ somando 0 e B_a ∈ [0, 1] nos valores de z
observados. `g` é o grupo (UF, região ou país) e `z_a` são covariáveis centradas da seção.
Estimação por mínimos quadrados ponderados pelo eleitorado, com a matriz de cada grupo
encolhida para a nacional; depois, IPF em cada seção para bater com os totais oficiais.

Efeito fixo de local de votação: o ajuste usa só as diferenças entre seções do mesmo local,

    y_a/E_a − ȳ_L ≈ B_aᵀ (x_a/E_a − x̄_L),

de modo que um efeito de contexto constante no local se cancela em vez de ser atribuído a um
grupo de origem. As diferenças só identificam B a menos de uma constante por coluna; o nível
de cada coluna vem da exigência de que a matriz reproduza o total do 2º turno do grupo.
"""
module InferenciaEcologica

using JuMP, Clarabel, LinearAlgebra, Random, Statistics

export REGIAO, estimar, sistema, fluxos, agrega, na_escala, Reamostrador, amostra, matriz, ipf

const REGIAO = Dict{String,String}()
for (reg, ufs) in ("Norte" => "ac am ap pa ro rr to", "Nordeste" => "al ba ce ma pb pe pi rn se",
                   "Centro-Oeste" => "df go ms mt", "Sudeste" => "es mg rj sp", "Sul" => "pr rs sc",
                   "Exterior" => "zz")
    for u in split(ufs)
        REGIAO[u] = reg
    end
end
const ORDEM_REGIOES = ["Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul", "Exterior"]

const K_GRUPO = 2.0e6   # encolhimento da matriz do grupo para a nacional, em eleitores
const N_PONTOS = 40     # pontos observados de z onde a matriz precisa ser válida
const IT_IPF = 60
const PISO_IPF = 1e-4   # probabilidade mínima de cada célula no ponto de partida do IPF


"Pontos de z (observados) em que exigimos 0 ≤ B ≤ 1: extremos de cada covariável e uma amostra."
function pontos_validade(Z::Matrix{Float64}, rng::AbstractRNG)
    idx = Set{Int}()
    for k in axes(Z, 2)
        push!(idx, argmin(Z[:, k]), argmax(Z[:, k]))
    end
    resto = setdiff(1:size(Z, 1), idx)
    union!(idx, resto[randperm(rng, length(resto))[1:min(length(resto), N_PONTOS - length(idx))]])
    Z[sort!(collect(idx)), :]
end

"Média (ponderada pelo eleitorado) de cada linha de A dentro do seu local (ids 1..L)."
function media_local(A::Matrix{Float64}, local_::Vector{Int}, n::Vector{Float64})
    L = maximum(local_)
    soma = zeros(L, size(A, 2))
    peso = zeros(L)
    @inbounds for i in eachindex(local_)
        l = local_[i]
        peso[l] += n[i]
        for c in axes(A, 2)
            soma[l, c] += A[i, c] * n[i]
        end
    end
    soma[local_, :] ./ peso[local_]
end

"Renumera ids para 1..L (preservando −1 como 'fora')."
function renumera(v::AbstractVector{<:Integer})
    mapa = Dict{Int,Int}()
    [x < 0 ? -1 : get!(mapa, x, length(mapa) + 1) for x in v]
end

"""
Matrizes normais G, H do ajuste (normalizadas pelo eleitorado total) e, com efeito fixo, o
vetor t e o alvo y_tot da restrição de totais (tᵀM = y_tot).

`local_`: id do local de cada seção, ou −1 para deixar a seção fora das diferenças. Locais
com uma só seção também ficam fora; todas as seções entram nos totais.
"""
function sistema(X::Matrix{Float64}, Y::Matrix{Float64}, Z, local_)
    n = vec(sum(X; dims=2))
    Y = Y .* (n ./ max.(vec(sum(Y; dims=2)), 1e-12))
    px, py = X ./ n, Y ./ n
    D = Z === nothing ? px : hcat(px, (px .* Z[:, c] for c in axes(Z, 2))...)
    if local_ === nothing
        w = n ./ sum(n)
        return D' * (D .* w), D' * (py .* w), nothing, nothing
    end
    ok = local_ .>= 0
    ids = renumera(local_[ok])
    cont = zeros(Int, maximum(ids; init=0))
    for i in ids
        cont[i] += 1
    end
    okk = findall(ok)[cont[ids] .> 1]
    lo = renumera(local_[okk])
    Dd = D[okk, :] .- media_local(D[okk, :], lo, n[okk])
    yd = py[okk, :] .- media_local(py[okk, :], lo, n[okk])
    w = n[okk] ./ sum(n)
    t = vec(sum(D .* n; dims=1)) ./ sum(n)
    Dd' * (Dd .* w), Dd' * (yd .* w), t, vec(sum(Y; dims=1)) ./ sum(n)
end

"""
Parâmetros empilhados ((1+K)·5 × 4): B⁽⁰⁾ em cima, depois um bloco por covariável.
`ref`/`k`: encolhimento quadrático para `ref` com peso k (em eleitores).
"""
function estimar(X::Matrix{Float64}, Y::Matrix{Float64}, Z=nothing; ref=nothing, k=0.0,
                 rng=MersenneTwister(0), local_=nothing, sis=nothing)
    n = vec(sum(X; dims=2))
    K = Z === nothing ? 0 : size(Z, 2)
    p = 5 * (1 + K)
    G, H, t, ytot = sis === nothing ? sistema(X, Y, Z, local_) : sis
    G = (G + G') / 2
    λ = ref === nothing ? 0.0 : k / sum(n)
    R = zeros(p, 4)
    ref === nothing || (R .= ref[1:p, :])
    resolve(otimizador) = begin
        m = Model(otimizador)
        set_silent(m)
        @variable(m, M[1:p, 1:4])
        @objective(m, Min, sum(0.5 * (M[:, j]' * G * M[:, j]) - H[:, j]' * M[:, j] for j in 1:4) +
                           0.5 * λ * sum((M[i, j] - R[i, j])^2 for i in 1:p, j in 1:4))
        for blk in 0:K, i in 1:5   # linhas somam 1 no bloco 0 e 0 nos blocos de covariável
            @constraint(m, sum(M[blk*5+i, :]) == (blk == 0 ? 1.0 : 0.0))
        end
        if t !== nothing           # efeito fixo: a matriz reproduz o total do 2º turno (3 colunas; a 4ª segue)
            for j in 1:3
                @constraint(m, dot(t, M[:, j]) == ytot[j])
            end
        end
        if K == 0
            @constraint(m, 0 .<= M .<= 1)
        else
            for z in eachrow(pontos_validade(Z, rng))
                B = [M[i, j] + sum(z[c] * M[5c+i, j] for c in 1:K) for i in 1:5, j in 1:4]
                @constraint(m, 0 .<= B .<= 1)
            end
            @constraint(m, 0 .<= M[1:5, :] .<= 1)
        end
        optimize!(m)
        is_solved_and_feasible(m) || error("solver: $(termination_status(m))")
        value.(M)
    end
    # Programação quadrática convexa (objetivo quadrático, restrições lineares), resolvida pelo
    # Clarabel, que é escrito em Julia e pode rodar em várias threads ao mesmo tempo.
    resolve(optimizer_with_attributes(Clarabel.Optimizer, "tol_gap_abs" => 1e-12, "tol_gap_rel" => 1e-12,
                                      "tol_feas" => 1e-12, "max_iter" => 500))
end

"Matriz de cada unidade (n×5×4), com B cortada em [0, 1] e linhas renormalizadas."
function matriz(θ::Matrix{Float64}, Z)
    n = Z === nothing ? 1 : size(Z, 1)
    M = repeat(reshape(θ[1:5, :], 1, 5, 4), n, 1, 1)
    if Z !== nothing
        for c in axes(Z, 2)
            M .+= Z[:, c] .* reshape(θ[5c+1:5c+5, :], 1, 5, 4)
        end
    end
    M = clamp.(M, 0, 1)
    M ./ max.(sum(M; dims=3), 1e-12)
end

"""
Fluxos por unidade (n×5×4) que batem com os totais da unidade nos dois turnos. Quando o
eleitorado registrado muda entre os turnos, o 2º turno é reescalado para o do 1º.
"""
function ipf(X::Matrix{Float64}, Y::Matrix{Float64}, M::Array{Float64,3}; it=IT_IPF, piso=PISO_IPF)
    n = size(X, 1)
    N = Array{Float64,3}(undef, n, 5, 4)
    Threads.@threads for a in 1:n   # cada seção é um problema 5×4 independente
        t = zeros(5, 4); x = zeros(5); y = zeros(4)
        @inbounds begin
            ex, ey = 0.0, 0.0
            for i in 1:5
                x[i] = X[a, i]; ex += x[i]
            end
            for j in 1:4
                ey += Y[a, j]
            end
            for j in 1:4
                y[j] = Y[a, j] * ex / max(ey, 1e-12)
            end
            for i in 1:5, j in 1:4
                t[i, j] = x[i] * (M[a, i, j] + piso)
            end
            for _ in 1:it
                for j in 1:4
                    c = 0.0
                    for i in 1:5
                        c += t[i, j]
                    end
                    f = y[j] / max(c, 1e-12)
                    for i in 1:5
                        t[i, j] *= f
                    end
                end
                for i in 1:5
                    r = 0.0
                    for j in 1:4
                        r += t[i, j]
                    end
                    f = x[i] / max(r, 1e-12)
                    for j in 1:4
                        t[i, j] *= f
                    end
                end
            end
            for i in 1:5, j in 1:4
                N[a, i, j] = t[i, j]
            end
        end
    end
    N
end

"""
Ajusta uma especificação e devolve os fluxos por unidade. `grupo`: :uf, :regiao ou nothing.
`local_`: ids de local de votação para o efeito fixo (ou nothing).
"""
function fluxos(X, Y, uf::Vector{String}, Z, grupo; local_=nothing, devolver_matrizes=false)
    sis = sistema(X, Y, Z, local_)
    nac = estimar(X, Y, Z; local_, sis)
    if grupo === nothing
        Ms = Z === nothing ? repeat(matriz(nac, nothing), size(X, 1), 1, 1) : matriz(nac, Z)
    else
        grupos = grupo === :uf ? uf : [REGIAO[u] for u in uf]
        # O encolhimento é dado em "eleitores equivalentes" do ajuste sem efeito fixo. Com efeito
        # fixo, cada eleitor carrega bem menos informação, então κ é reduzido na mesma proporção.
        k = K_GRUPO
        local_ === nothing || (k *= tr(sis[1]) / tr(sistema(X, Y, Z, nothing)[1]))
        Ms = Array{Float64,3}(undef, size(X, 1), 5, 4)
        for g in unique(grupos)
            s = findall(==(g), grupos)
            Zs = Z === nothing ? nothing : Z[s, :]
            th = estimar(X[s, :], Y[s, :], Zs; ref=nac, k, local_=local_ === nothing ? nothing : local_[s])
            Ms[s, :, :] = Zs === nothing ? repeat(matriz(th, nothing), length(s), 1, 1) : matriz(th, Zs)
        end
    end
    N = ipf(X, Y, Ms)
    devolver_matrizes ? (N, Ms) : N
end

"Soma os fluxos por recorte: Brasil, regiões e UFs."
function agrega(N::Array{Float64,3}, uf::Vector{String})
    soma(idx) = dropdims(sum(N[idx, :, :]; dims=1); dims=1)
    out = Dict{String,Matrix{Float64}}("Brasil" => soma(Colon()))
    regs = [REGIAO[u] for u in uf]
    for r in ORDEM_REGIOES
        out[r] = soma(regs .== r)
    end
    for u in unique(uf)
        out[uppercase(u)] = soma(uf .== u)
    end
    out
end

"""
Leva uma reamostragem ao eleitorado real de cada origem: o bootstrap por município sorteia
amostras de tamanho variável, mas o eleitorado de cada grupo de origem é conhecido; só as
taxas de transição são incertas.
"""
na_escala(boot, ponto) = Dict(nome => b ./ max.(sum(b; dims=2), 1) .* sum(ponto[nome]; dims=2) for (nome, b) in boot)

"Bootstrap por agrupamento: sorteia municípios (com reposição) dentro de cada UF."
struct Reamostrador
    rng::MersenneTwister
    por_uf::Vector{Vector{Vector{Int}}}
end
function Reamostrador(uf::Vector{String}, cluster::Vector{Int}, seed::Integer)
    por_uf = Vector{Vector{Int}}[]
    for u in sort(unique(uf))
        s = findall(==(u), uf)
        grupos = Dict{Int,Vector{Int}}()
        for i in s
            push!(get!(grupos, cluster[i], Int[]), i)
        end
        push!(por_uf, [grupos[c] for c in sort(collect(keys(grupos)))])
    end
    Reamostrador(MersenneTwister(seed), por_uf)
end
amostra(r::Reamostrador) = reduce(vcat, [g[rand(r.rng, 1:length(g))] for g in r.por_uf for _ in 1:length(g)])

end
