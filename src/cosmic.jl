const data_dir   = joinpath("..", "data")
const result_dir = joinpath("..", "result")
const model_dir  = joinpath("..", "model")

const test = false

@doc """ build machine learning datasets :: this seems parallel 
         ID_sample: [(Gene_name,feat),(),...]
         hash(Gene_name) --> col_idx 
""" ->
function build_ml_data()
    ml_dir = joinpath(data_dir, "ml")
    !isdir(ml_dir) && mkdir(ml_dir)
    
    ### genes:   all the genes then filter
    ### samples: all the samples then filter
    ### class: Primary_histology, Primary_site, sample_Source
    
    


end

@doc """ random build raw train, validation, evaluation dataset based on the raw_samples
""" ->
function split_samples(;ratios=(0.6,0.2,0.2))
    raw_samples = joinpath(data_dir, "raw_samples")
    sample_fls = try
        readdir(raw_samples)
    catch 
        @warn("$raw_samples is not exits, program terminated") #TODO
        exit(-1)
    end
    num_samples = length(sample_fls) - 1 
    @assert length(readdir(raw_samples)) > 0
    
    tr_val_point  = round(Int64, num_samples * ratios[1])
    val_te_point  = round(Int64, num_samples * (ratios[1] + ratios[2]))
    
    shuffle!(sample_fls)
    train_samples = sample_fls[1:tr_val_point]
    val_samples   = sample_fls[tr_val_point+1:val_te_point]
    te_samples    = sample_fls[val_te_point+1:end]

    sample_meta   = joinpath(data_dir, "sample_meta")
    if !isdir(sample_meta)
        mkdir(sample_meta)
    end
    train_fl,val_fl,test_fl = "train.csv","val.csv","test.csv"
    for (fl,samples) in [(train_fl,train_samples), (val_fl,val_samples),(test_fl, te_samples)]
        fl = joinpath(sample_meta, fl)
        writecsv(fl, samples)
    end
end

@doc """ build raw samples to disk: ../data/raw_samples from raw cosmic file
""" ->
function build_raw_samples(cosmic_path::ASCIIString)
    raw_samples = joinpath(data_dir,"raw_samples")
    if !isdir(raw_samples)
        mkdir(raw_samples)
    else
        println("removing files under $raw_samples")
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


