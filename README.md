# pq_inferencia_ecologica

Inferência ecológica do fluxo de votos entre o 1º e o 2º turno da eleição presidencial de 2026, por seção eleitoral, em Julia. Código de apoio ao artigo do [Pílulas de Quant](https://pilulasdequant.substack.com) sobre inferência ecológica com dados agregados e suas aplicações em finanças quantitativas.

O voto é secreto e o TSE publica apenas totais por seção. O modelo estima a matriz de transição (de cada 100 eleitores de cada grupo do 1º turno, quantos foram para cada grupo do 2º) a partir dos resultados das cerca de 500 mil seções, com faixas de incerteza calibradas por simulação.

## Modelo

- **Regressão ecológica com restrições** (Goodman): mínimos quadrados ponderados pelo eleitorado, com taxas entre 0 e 1 e linhas somando 1, resolvidos como programação quadrática (JuMP + HiGHS).
- **Covariáveis**: a matriz varia, dentro de cada estado, com a escolaridade e a idade do eleitorado da seção (perfil do eleitorado do TSE).
- **Encolhimento**: a matriz de cada estado é puxada para a nacional (Bayes empírico).
- **Efeito fixo de local de votação**: o ajuste usa só as diferenças entre seções da mesma escola, de modo que um efeito de contexto do bairro não é atribuído a um grupo de eleitores. O nível de cada coluna vem da exigência de que a matriz reproduza o resultado do 2º turno em cada estado.
- **Ajuste proporcional iterativo (IPF)** em cada seção, para bater com os totais oficiais.
- **Incerteza**: bootstrap por município, cinco versões do modelo e calibração da faixa em eleições simuladas com verdade conhecida.
- **Saldos**: ida menos volta para cada par de grupos, com faixa própria. A ida e a volta separadas são pouco identificadas; o saldo é bem mais firme.

A formulação completa está nos comentários de `src/InferenciaEcologica.jl` e de cada script.

## Como reproduzir

Julia 1.10 ou mais nova.

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/baixar.jl 2026           # dados abertos do TSE em data/brutos/2026 (~1,5 GB)
julia --project=. scripts/baixar.jl 2022           # base das eleições simuladas da validação
julia -t 8 --project=. scripts/preparar.jl 2026    # data/secoes/2026.arrow, uma linha por seção
julia -t 8 --project=. scripts/preparar.jl 2022
julia -t 8 --project=. scripts/validar.jl 24 30    # resultados/validacao.json (calibração das faixas)
julia -t 8 --project=. scripts/fluxos.jl 2026 400  # resultados/fluxos-secao-2026.json
```

Rode a validação antes dos fluxos: a faixa publicada usa a calibração gravada em `resultados/validacao.json` (sem o arquivo, sai a faixa bruta). Para usar arquivos do TSE já baixados em outra pasta, defina `TSE_BRUTOS` antes de `preparar.jl`.

`scripts/compara.jl` compara dois arquivos de fluxos no mesmo formato, por exemplo este resultado e o da implementação original em Python.

## Resultados

Brasil, 2026, em milhões de eleitores, com a faixa calibrada (`resultados/fluxos-secao-2026.json`):

| 1º turno → 2º turno | Estimativa | Faixa |
|---|---|---|
| Outros candidatos → Lula | 4,01 | 3,48 a 4,18 |
| Outros candidatos → Flávio | 3,63 | 2,94 a 3,77 |
| Abstenção → Lula | 4,01 | 3,51 a 4,31 |
| Abstenção → Flávio | 1,78 | 0,51 a 2,20 |
| Brancos e nulos → Lula | 2,68 | 2,55 a 3,30 |
| Brancos e nulos → Flávio | 1,17 | 0,69 a 1,35 |
| Flávio → Lula | 2,06 | 0,00 a 2,27 |
| Lula → Flávio | 0,61 | 0,00 a 2,13 |

A validação (`resultados/validacao.json`) usa 96 eleições simuladas sobre 146 mil seções de 2022. A faixa calibrada (c = 1,6 e d = 1,4%) contém a verdade em 96% dos fluxos relevantes sem efeito de contexto (cenário A), 92% com heterogeneidade não observada (B) e 81% com contexto por local de votação (D). Com contexto que varia de seção para seção dentro da mesma escola (C), nenhuma versão acerta a linha de brancos e nulos.

## Conferência com a implementação em Python

O painel da apuração de 2026 usa uma implementação anterior em Python. As duas foram comparadas:

- **Tabela por seção:** idêntica nas 499.207 seções de 2026 e nas 472.028 de 2022.
- **Estimativas:** diferem em até 57 mil eleitores no total nacional. Sem covariáveis, a diferença cai para cerca de 10 mil e vem do otimizador do Python (SLSQP), que para antes do ótimo (o Clarabel chega a um valor menor da função objetivo). Na versão principal soma-se o sorteio dos pontos em que a restrição 0 ≤ B ≤ 1 é imposta, que só ele move a matriz nacional em até 25 mil eleitores.
- **Validação:** cobertura calibrada de 96%, 92%, 44% e 81% nos cenários A a D, contra 97%, 94%, 42% e 82% no Python. As diferenças vêm dos geradores aleatórios.
- **Tempo, 8 threads:** preparo dos dados em cerca de 50 s, fluxos com 400 reamostragens em 5 min (Python: 18 min) e validação em 13 min (Python: 84 min).

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `src/InferenciaEcologica.jl` | núcleo do modelo: estimação com restrições, efeito fixo de local, IPF, bootstrap por município |
| `scripts/baixar.jl` | download dos arquivos por seção do TSE |
| `scripts/preparar.jl` | tabela por seção: votos nos dois turnos e perfil do eleitorado |
| `scripts/fluxos.jl` | versões do modelo, bootstrap, faixas calibradas e saldos |
| `scripts/validar.jl` | eleições simuladas em quatro cenários e calibração das faixas |
| `scripts/compara.jl` | comparação entre duas saídas |
| `resultados/` | saídas JSON |

Os dados brutos e a tabela por seção ficam fora do git (são grandes e públicos no portal de dados abertos do TSE).
