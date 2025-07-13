"""
    whiten(data::AbstractMatrix) -> Tuple{Matrix, Matrix}

Performs whitening on the given data matrix.

Whitening transforms the input signals so that each signal has unit variance and all signals are uncorrelated.
This is a common preprocessing step before applying Independent Component Analysis (ICA).

# Arguments
- `data`: An abstract matrix of size `(n_samples, n_signals+1)` where the first column is time and the remaining `n_signals` columns are the observed signals.

# Returns
- `data_white::Matrix`: A matrix of the same size as `data`, where the first column is the original time vector and the remaining columns are the whitened signals.
- `W_white::Matrix`: The whitening transform matrix of size `(n_signals, n_signals)` such that `X_white = X_centered * W_white`.

"""
function whiten(data::AbstractMatrix)
    time = @view data[:, 1]
    X = @view data[:, 2:end]

    X_mean = mean(X, dims=1)
    X_centered = X .- X_mean

    U, S, Vt = svd(X_centered, full=true)

    W_white = Vt' * Diagonal(1 ./ S)

    X_white = X_centered * W_white

    data_white = hcat(time, X_white)

    return data_white, W_white
end