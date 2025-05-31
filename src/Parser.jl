using DelimitedFiles
using Plots

"""
    ReadDatasetFromDatFile(Path :: String,tabdelimited :: Bool = true)

Reads .dat File from a File Location. 
returns 9x2497 Matrix {Float64} with oure dataset
"""
function ReadDatasetFromDatFile(Path :: String)
    return readdlm(Path,Float64)
end



