"""
    read_dataset_from_dat(Path :: String)

Reads .dat File from a File Location. 
returns 9x2497 Matrix {Float64} with our dataset
"""
function read_dataset_from_dat(Path :: String)
    return readdlm(Path,Float64)
end



