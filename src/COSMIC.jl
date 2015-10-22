module COSMIC

using Logging
@Logging.configure(level=DEBUG)
using DataFrames
using Compat
#using Docile
using HDF5

export 
    build_raw_samples



include("utils.jl")
include("cosmic.jl")



function __init__()
#     @pline "1"
end

end # module

