module ICAforECGrecordings

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
    centered = data .- mean(data, dims=1)
    cov_matrix = cov(centered, dims=1)
    eigvals, eigvecs = eigen(cov_matrix)
    whitening_matrix = inv(sqrt.(Diagonal(eigvals))) * eigvecs'
    whitened = centered * whitening_matrix'
    return whitened
end
end
