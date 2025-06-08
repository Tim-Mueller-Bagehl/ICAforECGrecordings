module ICAforECGrecordings
using LinearAlgebra
using Statistics

# Write your package code here.
export plot_dataset
include("functions.jl")

export ReadDatasetFromDatFile
include("Parser.jl")

export whiten

"""
whiten(data::Matrix{Float64}) -> Matrix{Float64}

Performs whitening on the given data matrix.

Whitening transforms the input signals so that each signal has unit variance and all signals are uncorrelated.
This is a common preprocessing step before applying Independent Component Analysis (ICA).

# Arguments
- `data`: A matrix of size (n_samples, n_signals) where each column is a signal.

# Returns
- A whitened data matrix of the same size, where the signals are zero-mean, uncorrelated, and have unit variance.

"""
function whiten(data::Matrix{Float64})
    time = data[:, 1]
    X = data[:, 2:end]

    X_mean = mean(X, dims=1)
    X_centered = X .- X_mean

    U, S, Vt = svd(X_centered; full=true)

    W_white = Vt' * Diagonal(1.0 ./ S)

    X_white = X_centered * W_white

    data_white = hcat(time, X_white)

    return data_white
end
end