const data_dir   = "../data"
const result_dir = "../result"
const model_dir  = "../model"

const test = false



@doc """ build raw train, validation, evaluation dataset based on the raw_samples
""" ->
function build_datasets()


end
#@doc """ build raw samples to disk: ../data/raw_samples from raw cosmic file
#""" ->
function build_raw_samples(cosmic_path::ASCIIString)
    raw_samples = joinpath(data_dir,"raw_samples")
    if !isdir(raw_samples)
        mkdir(raw_samples)
    else
        @info("removing files under $raw_samples")
        rm(raw_samples,recursive=true)
        mkdir(raw_samples)
    end
    f = open(cosmic_path)
    header = readline(f)
    writedlm(joinpath(raw_samples, "header.tsv"), header)
    for line in eachline(f)
        row = split(line,'\t')
        sample_name = strip(row[6]) # #row[6] is ID_sample
        f_sp = open(joinpath(raw_samples,string(sample_name,".tsv")), "a+")
        write(f_sp, row)
        close(f_sp)
    end
    close(f)
end

function build_dir()
    for dir in (data_dir,result_dir,model_dir)
        if !isdir(dir)
            println("buiding directory $data_dir")
            mkdir(dir)
        end
    end
end


