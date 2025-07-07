module ICAforECGrecordings
using LinearAlgebra: tranpose, sqrt, I, norm
using DelimitedFiles: readdlm
using Statistics: norm, mean
using Plots: plot, plot!, xlabel!


# Write your package code here.
export whiten
include("Preprocessing.jl")

export plot_dataset
include("Visualization.jl")

export ReadDatasetFromDatFile
include("Parser.jl")


export shibbs
include("Shibbs.jl")

export jade
include("Jade.jl")

end