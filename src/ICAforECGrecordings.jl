module ICAforECGrecordings

# Write your package code here.
export plot_dataset
include("functions.jl")

export ReadDatasetFromDatFile
include("Parser.jl")

end
