# Compara o resultado desta implementação com outro arquivo de fluxos no mesmo formato
# (por exemplo, o gerado pela implementação original em Python do painel da apuração).
#
# Uso: julia --project=. scripts/compara.jl resultados/fluxos-secao-2026.json outro.json

using JSON3, Printf

a, b = (JSON3.read(read(f, String)) for f in ARGS[1:2])
mat(x) = [Float64(x[i][j]) for i in 1:5, j in 1:4]
println("Recorte        maior |Δ| estimativa   maior |Δ| ponta da faixa   (milhares de eleitores)")
pior = (0.0, "")
for nome in vcat(["Brasil", "Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul"], sort([String(k) for k in keys(a.niveis) if length(String(k)) == 2]))
    na, nb = a.niveis[Symbol(nome)], b.niveis[Symbol(nome)]
    de = maximum(abs.(mat(na.fluxos) .- mat(nb.fluxos)))
    df = max(maximum(abs.(mat(na.faixa[1]) .- mat(nb.faixa[1]))), maximum(abs.(mat(na.faixa[2]) .- mat(nb.faixa[2]))))
    de > pior[1] && (global pior = (de, nome))
    nome in ("Brasil", "Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul") && @printf("%-14s %10.1f %26.1f\n", nome, de / 1e3, df / 1e3)
end
@printf("\nMaior diferença de estimativa entre todos os recortes: %.1f mil eleitores (%s)\n", pior[1] / 1e3, pior[2])
println("\nBrasil, estimativa (milhões): esta implementação | outra")
O = ["Lula", "Outros", "Adversário", "Brancos/nulos", "Abstenção"]
for i in 1:5
    @printf("  %-14s %s\n", O[i], join([@sprintf("%6.3f | %6.3f", a.niveis.Brasil.fluxos[i][j] / 1e6, b.niveis.Brasil.fluxos[i][j] / 1e6) for j in 1:4], "   "))
end
if haskey(a.niveis.Brasil, :saldos) && haskey(b.niveis.Brasil, :saldos)
    println("\nSaldos, Brasil (milhões): saldo [faixa]")
    for (k, s) in pairs(a.niveis.Brasil.saldos)
        t = b.niveis.Brasil.saldos[k]
        @printf("  %s  %6.2f [%5.2f, %5.2f]  |  %6.2f [%5.2f, %5.2f]\n", k, s.saldo / 1e6, s.faixa[1] / 1e6, s.faixa[2] / 1e6, t.saldo / 1e6, t.faixa[1] / 1e6, t.faixa[2] / 1e6)
    end
end
