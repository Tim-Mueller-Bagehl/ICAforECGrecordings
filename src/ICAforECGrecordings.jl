module ICAforECGrecordings
using LinearAlgebra: transpose, sqrt, I, norm, svd, Diagonal, pinv ,diag, diagm, qr, rank
using DelimitedFiles: readdlm
using Statistics: norm, mean
using Plots: plot, plot!, xlabel!
using Random: randn



# Write your package code here.
include("Preprocessing.jl")
include("Visualization.jl")
include("Parser.jl")

# Algos
include("Shibbs.jl")
include("Jade.jl")
include("Picard.jl")


@doc """
    load_example_data()
Load example data from a .dat file for testing purposes.
Returns a matrix containing the data.
"""
function load_example_data()
    datpath = joinpath(@__DIR__, "..", "data", "foetal_ecg.dat")
    data = read_dataset_from_dat(datpath)
    return data
end

export whiten, plot_dataset, read_dataset_from_dat, shibbs, jade, load_example_data, picard


end 