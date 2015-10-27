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
function col(atom::Atom)
    atom.idx_col
end
function field_val(atom::Atom)
    atom.field_val
end
function convert(Sample, atom::Atom)
    Sample(atom.idx_row, Int64[atom.idx_col], ASCIIString[field_val(atom)])
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
function field_val(sample::Sample)
    sample.field_val
end
function row(sample::Sample)
    sample.idx_row
end
function col(sample::Sample)
    sample.idx_col
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
    Sample(row(atom_1), col(atom_1,atom_2),field_val(atom_1,atom_2))
end
function Sample(sample::Sample, atom::Atom)
    @assert row(sample) == row(atom)
    idx_cols = vcat(col(sample),col(atom))
    field_vals = vcat(field_val(sample),field_val(atom))
    Sample(row(sample),idx_cols,field_vals)
end
function Sample(sample_1::Sample, sample_2::Sample)
    @assert row(sample_1) == row(sample_2)
    idx_cols = vcat(col(sample_1),col(sample_2))
    field_vals = vcat(field_val(sample_1),field_val(sample_2))
    Sample(row(sample_1),idx_cols,field_vals)
end

function show(sample::Sample)
    println(sample.idx_row,sample.idx_col)
    println(sample.field_val)
end
function tet()
end
function row(samples::Array{Sample})
    map(row, samples)
end
type OrderedSampleMatrix
    num_samples
    samples::Array{Sample} # sorted idx_IDsample
end
function OrderedSampleMatrix(samples::Array{Sample}) # for null ordered sample matrix
    idxs = row(samples)
    num_samples = length(samples)
    perm = sortperm(idxs)
    OrderedSampleMatrix(num_samples, samples[perm])
end
function start(orspmat::OrderedSampleMatrix)
    1
end
function start(sp::Sample)
    1
end
function next(sp::Sample, state)
    sp.samples[state], state+1
end
function done(orspmat::OrderedSampleMatrix, state)
    state > orspmat.num_samples
end
function convert(::Type{Array{Sample,1}}, samples::Array{Any,1})
    ret = Sample[]
    for sample in samples
        push!(ret, sample)
    end
    ret
end

#@doc """  build selected fields into OrderedSampleMatrices at field_dir
#""" ->
function build_field_datasets()
    header = readdlm(joinpath(sample_meta, "header.tsv"), '\t', ASCIIString)
    header = header[:,1]
    @show header
    deleteat!(header, 1)
    field_names = [(name,idx) for (idx,name) in enumerate(header)]
    map(build_field_dataset, field_names)
end
function build_field_dataset(field_meta::Tuple{ASCIIString,Int64})
    build_field_dataset(field_meta[1],field_meta[2])
end
#@doc """ the simplist way of build machine learning dataset only using one feature: Mutation_ID
#         ID_sample: [(Gene_name,feat),(),...]
#         hash(Gene_name) --> col_idx 
#         hash(Sample_ID) --> row_idx
#""" ->
function build_field_dataset(field_name::ASCIIString, idx_field::Int64)
    field_dir = joinpath(data_dir, "field")
    !isdir(field_dir) && mkdir(field_dir)
    
    # create hashes for sample_ID and single_field
    
    ### genes:   all the genes then filter
    ### samples: all the samples then filter
    ### field: single_field
    sample_meta   = joinpath(data_dir, "sample_meta")
    !isdir(sample_meta) && @error("can't find sample_meta")

    field_vals = readcsv(joinpath(sample_meta, string(field_name, ".csv")))
    gene_names = readcsv(joinpath(sample_meta, "Gene name.csv"), ASCIIString)
    datasets = ["train.csv","val.csv","test.csv"]
    num = 0
    for dataset in datasets
        println("dataset is ", dataset)
        sample_fls = readcsv(joinpath(sample_meta,dataset))
        ID_samples = map(x->x[1:end-4],sample_fls)
        hash_sample = Dict()
        for idx = 1:length(sample_fls)
            hash_sample[ID_samples[idx]] = idx
        end

        #the name of columns is Gene_name

        hash_gene = Dict() 
        for idx = 1:length(gene_names)
            hash_gene[gene_names[idx]] = idx
        end
        @show hash_gene["39340"]
        function merge_atom(atom_1::Atom,atom_2::Atom)
            Sample(atom_1, atom_2) 
        end
        function merge_atom(atom_1::Sample,atom_2::Atom)
            Sample(atom_1, atom_2) 
        end
        function merge_atom(atom_1::Sample,atom_2::Sample)
            Sample(atom_1, atom_2) 
        end

        samples = @parallel vcat for sample_fl in sample_fls
            f = open(joinpath(raw_samples, sample_fl))
            lines = readlines(f)
            function atom_mapper(line)
                line = strip(line,'\n')
                row = split(line, '\t')
                @assert length(row) == 34
                # row[6] is ID_sample
                row_idx = try 
                            hash_sample[row[6]]
                          catch e
                            println("==========Error Info===========")
                            println(e)
                            println(row)
                          end
                        
                col_idx = try 
                            hash_gene[strip(row[1])]
                          catch e
                            println("==========Error Info===========")
                            println(e)
                            println(row)
                          end
                           
                # Mutation_ID:#17 row_idx col_idx
                atom = Atom(row_idx, col_idx, row[idx_field])
            end
            sample = reduce(merge_atom, pmap(atom_mapper, lines)) # short filea
            if typeof(sample) <: Atom
                sample = convert(Sample, sample)
            end
            close(f)
            sample
        end
        @debug("here ========")
        sample_matrix = OrderedSampleMatrix(samples)
        dataset_dir = joinpath(field_dir, filename(dataset))
        isdir(dataset_dir) || mkdir(dataset_dir)
        dataset_name = joinpath(dataset_dir, string(field_name,".tsv"))
        @debug("bugs above ====")
        jldopen(dataset_name, "w") do file
            println("save $dataset into file $dataset_name")
            write(file, "sample_matrix", sample_matrix)
        end
    end
end

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
    write(joinpath(sample_meta, "header.tsv"), split(strip(header, '\t')))
    for line in eachline(f)
        row = split(line,'\t')
        sample_name = strip(row[6]) # #row[6] is ID_sample
        f_sp = open(joinpath(raw_samples,string(sample_name,".tsv")), "a+")
        write(f_sp, line)
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
