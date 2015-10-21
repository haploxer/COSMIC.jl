module COSMIC

using Logging
@Logging.configure(level=DEBUG)
using DataFrames
using Compat
using Docile


export 
    test,
    build_raw_data



include("cosmic.jl")
include("utils.jl")



function __init__()
#    build_dir()
end

end # module

