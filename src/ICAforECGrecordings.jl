module ICAforECGrecordings
using LinearAlgebra: transpose, sqrt, I, norm, svd, Diagonal, pinv
using DelimitedFiles: readdlm
using Statistics: norm, mean
using Plots: plot, plot!, xlabel!


# Write your package code here.
export whiten
include("Preprocessing.jl")

export plot_dataset
include("Visualization.jl")

export read_dataset_from_dat
include("Parser.jl")


export shibbs
include("Shibbs.jl")

export jade
include("Jade.jl")


end