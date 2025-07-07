"""
    whiten(data::AbstractMatrix{Float64}) -> AbstractMatrix{Float64}

Performs whitening on the given data matrix.

Whitening transforms the input signals so that each signal has unit variance and all signals are uncorrelated.
This is a common preprocessing step before applying Independent Component Analysis (ICA).

# Arguments
- `data`: A matrix of size (n_samples, n_signals) where each column is a signal.

# Returns
- A matrix of the same size as `data`, where the first column is time and the remaining columns are the whitened signals.

"""
function whiten(data::AbstractMatrix{Float64})
    time = data[:, 1]
    X = data[:, 2:end]

    X_mean = mean(X, dims=1)
    X_centered = X .- X_mean

    U, S, Vt = svd(X_centered; full=true)

    W_white = Vt' * Diagonal(1 ./ S)

    X_white = X_centered * W_white

    data_white = hcat(time, X_white)

    return data_white
end