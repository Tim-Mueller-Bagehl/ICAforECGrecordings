module ICAforECGrecordings
using LinearAlgebra: I, eigen, Symmetric, dot, inv, Diagonal
using Statistics: mean

# Write your package code here.
export whiten
include("Preprocessing.jl")

export plot_dataset
include("Visualization.jl")


export ReadDatasetFromDatFile
include("Parser.jl")

export jade
include("Jade.jl")

end