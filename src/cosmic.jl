const data_dir    = joinpath("..", "data")
const result_dir  = joinpath("..", "result")
const model_dir   = joinpath("..", "model")
const raw_samples = joinpath(data_dir, "raw_samples")
const sample_meta = joinpath(data_dir, "sample_meta")
const raw_field_datasets = joinpath(data_dir,"raw_field_datasets")

const miss = "MISSING"
const test = false

immutable Atom
    idx_row::Int64
    idx_col::Int64
    field_val::ASCIIString
end
function row(atom::Atom)
    atom.idx_row
end
function col(atom_1::Atom,atom_2::Atom)
    Int64[atom_1.idx_col, atom_2.idx_col]
end
function field_val(atom_1::Atom, atom_2::Atom)
    ASCIIString[atom_1.field_val, atom_2.field_val]
end
function Atom(idx_row,idx_col,field_val::ASCIIString)
    @assert idx_row > 0
    @assert idx_col > 0
    if length(field_val) == 0 
        field_val = miss
    end
    Atom(idx_row, idx_col, field_val)
end

immutable Sample
    idx_row::Int64
    idx_col::Array{Int64,1}
    field_val::Array{ASCIIString,1}
end
function Sample(idx_row,idx_col,field_val)
    @assert idx_row > 0
    @assert idx_col > 0
    if length(field_val) == 0  
        field_val = miss
    end

    Sample(idx_row, idx_col, field_val)
end
function Sample(atom_1::Atom, atom_2::Atom)
    @assert row(atom_1) == row(atom_2)
    Sample(idx_row,col(atom_1,atom_2),field_val(atom_1,atom_2))
end
function row(sample::Sample)
    sample.idx_row
end
function row(samples::Array{Sample})
    map(row, samples)
end
immutable OrderedSampleMatrix
    samples::Array{Sample} # sorted idx_IDsample
end
function OrderedSampleMatrix(samples::Array{Sample}) # for null ordered sample matrix
    idxs = row(samples)
    perm = sortperm(idxs)
    OrderedSampleMatrix(samples[perm])
end

@doc """  build selected fields into OrderedSampleMatrices at field_dir
""" ->
@acc function build_field_datasets()
    const field_names = [("Mutation_ID",17),("FATHMM_score",34)]
    map(build_field_dataset, field_names)
end
function build_field_dataset(field_meta::Tuple{ASCIIString,Int64})
    build_field_dataset(field_meta[1],field_meta[2])
end
@doc """ the simplist way of build machine learning dataset only using one feature: Mutation_ID
         ID_sample: [(Gene_name,feat),(),...]
         hash(Gene_name) --> col_idx 
         hash(Sample_ID) --> row_idx
""" ->
@acc function build_field_dataset(field_name::ASCIIString, idx_field::Int64)
    field_dir = joinpath(data_dir, "field")
    !isdir(field_dir) && mkdir(field_dir)
    
    # create hashes for sample_ID and single_field
    
    ### genes:   all the genes then filter
    ### samples: all the samples then filter
    ### field: single_field
    sample_meta   = joinpath(data_dir, "sample_meta")
    !isdir(sample_meta) && @error("can't find sample_meta")

    field_vals = readcsv(joinpath(sample_meta, field_name, ".csv"))
    for dataset in readdir(sample_meta)

        sample_fls = readcsv(joinpath(sample_meta,dataset))
        hash_sample = Dict()
        for idx = 1:length(sample_fls)
            hash_sample[sample_fls[idx]] = idx
        end

        hash_field = Dict()
        for idx = 1:length(field_vals)
            hash_field[field_vals[idx]] = idx
        end

        function merge_atom(atom_1,atom_2)
            Sample(atom_1, atom_2) 
        end

        samples = @parallel for sample_fl in sample_fls
            f = open(joinpath(raw_samples, sample_fl))

            function atom_mapper(line)
                row = split(line, '\t')
                # row[6] is ID_sample
                row_idx = hash_sample(row[6])
                col_idx = hash_field(row[idx_field])
                # Mutation_ID:#17 row_idx col_idx
                Atom(row_idx, col_idx, row[idx_field])
            end
            sample = reduce(merge_atom, pmap(atom_mapper, readlines(f))) # short filea

            close(f)
            sample
        end
        sample_matrix = OrderedSampleMatrix(samples)
        dataset_name = joinpath(field_dir, filename(dataset), field_name)
        writedlm(dataset_name, sample_matrix)
    end

end

#= it seems that his function is deprecated
@doc """ build machine learning datasets :: this seems parallel 
         ID_sample: [(Gene_name,feat),(),...]
         hash(Gene_name) --> col_idx 
""" ->
function build_ml_data(;sample_fields=("Gene_name","Gene_CDS_length","Primary_site",
                        "Primary_histology","Histology_subtype_1","Mutation_ID","Mutation_Description"
                        ,"Mutation_zygosity","LOH","SNP","Sample_source","Tumour_origin","Age"))
    ml_dir = joinpath(data_dir, "ml")
    !isdir(ml_dir) && mkdir(ml_dir)
    
    ### genes:   all the genes then filter
    ### samples: all the samples then filter
    ### class: Primary_histology, Primary_site, sample_Source
    sample_meta   = joinpath(data_dir, "sample_meta")
    !isdir(sample_meta) && @error("can't find sample_meta")
    for dataset in readdir(sample_meta)
        for sample_fl in readcsv(joinpath(sample_meta,dataset)) # @parallel
            f = open(joinpath(raw_samples,sample_fl))
            for line in eachline(f)
                row = split(line, '\t')
                @warn("Implement of this function not finished yet")
            end
            close(f)
        end
    end

end
=#


@doc """ Since cosmic raw data file is large, we tend to build all the fields with only one exaust search
""" ->
function build_fields_meta(cosmic_path::ASCIIString)
    f = open(cosmic_path)
    header = readline(f)
    header = split(header, '\t')
    data   = readdlm(f, ASCIIString)
    map!(x -> length(x) == 0 ? miss : x, data)
    
    if isdir(sample_meta) 
        
    else
        mkdir(sample_meta)
    end
    for (idx,field) in enumerate(header) #@simd
        writecsv(joinpath(sample_meta, string(strip(header[idx]),".csv")), sort(unique(data[:,idx])))
    end

    #=
    fields = header[]
    fl_meta = joinpath(sample_meta, name_meta, ".csv")
#   meta_io = open(meta_fl, "w")
    ret = Array{ASCIIString}[]
    for line in eachline(f) # can't parallel since data is too large
        row = split(line, '\t')
        # Gene_name:#1, Gene_CDS_length:#3
        push!(ret, ASCIIString[row[i] for i in idx_meta]')
    end
    ret = reduce(vcat, unique(ret))
    writecsv(fl_meta, ret)
    =#
end

#function build_field()
#end

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

function filename(file::ASCIIString)
    split(file,".")[1]
end
