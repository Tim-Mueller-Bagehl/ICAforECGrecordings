"""
    read_dataset_from_dat(Path :: String)
Reads a dataset from a .dat file.
# Arguments
- `Path::String`: The path to the .dat file.
# Returns
- A matrix of type `Float64` containing the data read from the file.
"""
function read_dataset_from_dat(Path :: String)
    return readdlm(Path,Float64)
end



