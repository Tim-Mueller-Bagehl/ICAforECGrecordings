using DelimitedFiles
"""
    ReadDatasetFromDatFile(Path :: String)

Reads .dat File from a File Location. 
returns 9x2497 Matrix {Float64} with our dataset
"""
function ReadDatasetFromDatFile(Path :: String)
    return readdlm(Path,Float64)
end



