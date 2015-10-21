if "../src" ∉  LOAD_PATH
    push!(LOAD_PATH,"../src")
end

import COSMIC

#const data_dir = "/home/guo/haplox/HapBrain/Tuhao/data"
const cosmic = "/home/guo/haplox/HapBrain/data/cosmic/GRch37/v74/CosmicCompleteExport.tsv"
const test_cosmic = "/home/guo/haplox/HapBrain/data/cosmic/10k.tsv"
#const result_dir = "/home/guo/haplox/HapBrain/Tuhao/result"
#const baseinfo = joinpath(result_dir,"baseinfo.txt")
const gene_census = "/home/guo/haplox/HapBrain/data/cosmic/cancer_gene_census.csv"
#const gene_tuhao = joinpath(data_dir, "tuhao_genes.csv")

#const debug = false
#const test = false


COSMIC.build_raw_samples(COSMIC.@test cosmic)

