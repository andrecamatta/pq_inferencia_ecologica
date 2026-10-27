# Baixa do portal de dados abertos do TSE os arquivos por seção de um ano.
#
# Salva em data/brutos/<ano>/ (fora do git). O CDN do TSE às vezes entrega a versão em cache
# mesmo depois de o arquivo mudar: o script pede sempre a versão atual e guarda a data dela ao
# lado do arquivo; se o arquivo mudou no TSE, baixa de novo.
#
# Uso: julia --project=. scripts/baixar.jl 2026

using Downloads

const BASE = "https://cdn.tse.jus.br/estatistica/sead/odsele"
const UFS = split("AC AL AM AP BA CE DF ES GO MA MG MS MT PA PB PE PI PR RJ RN RO RR RS SC SE SP TO ZZ")

function baixar(url, destino)
    url *= "?v=$(round(Int, time()))"
    sem_cache = ["Cache-Control" => "no-cache"]
    cab = Dict(lowercase(k) => v for (k, v) in Downloads.request(url; method="HEAD", headers=sem_cache).headers)
    versao, tam = get(cab, "last-modified", ""), parse(Int, cab["content-length"])
    marca = destino * ".versao"
    if isfile(destino) && filesize(destino) == tam && isfile(marca) && read(marca, String) == versao
        return
    end
    Downloads.download(url, destino; headers=sem_cache)
    write(marca, versao)
    println("ok $(basename(destino)) ($(round(Int, tam / 1e6)) MB)")
end

function main(ano)
    pasta = joinpath(dirname(@__DIR__), "data", "brutos", ano)
    mkpath(pasta)
    arquivos = ["votacao_secao/votacao_secao_$(ano)_BR.zip", "detalhe_votacao_secao/detalhe_votacao_secao_$(ano).zip"]
    append!(arquivos, ["perfil_eleitor_secao/perfil_eleitor_secao_$(ano)_$(u).zip" for u in UFS])
    for a in arquivos
        baixar("$BASE/$a", joinpath(pasta, basename(a)))
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__()) && main(ARGS[1])
