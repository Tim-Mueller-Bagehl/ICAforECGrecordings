"""
    jade(data::AbstractMatrix{<:Real}) -> Matrix{Float64}

Perform Blind Source Separation on multichannel time-series data using the JADE (Joint
Approximate Diagonalization of Eigen-matrices) algorithm.

# Arguments
- `data_w::AbstractMatrix{<:Real}`  
  A matrix of size `(n_samples, n_signals + 1)` where the first column is time in seconds
  and the remaining columns are the observed signals.
- `W_white::AbstractMatrix{<:Real}`
  The whitening matrix of size `(n_signals, n_signals)` used to preprocess the data.

# Returns
- A matrix of size `(n_samples, n_signals + 1)` where the first column is time and
  the remaining columns are the estimated source signals after applying JADE.
"""
function jade(data_w::AbstractMatrix, W::AbstractMatrix) 
    
    time = data_w[:, 1]
    Xs = data_w[:, 2:end]
    X = Xs'
    _,T = size(X)

    CM = fourth_order_cumulant_matrices(Matrix(X))
    V = joint_diagonalization_new(CM, T, 2)
    
    B = extract_sources(V, W)
    #B = V' * W
    S = B * X
    iS = S'
    return hcat(time, iS[:,1:2])
end