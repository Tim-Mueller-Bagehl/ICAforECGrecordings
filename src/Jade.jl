"""
    separate_sources(V::AbstractMatrix, W::AbstractMatrix) -> B::Matrix{Float64}

Compute a scaled and permuted separation matrix `B` given a demixing matrix `V` and
a mixing matrix `W`, ordering output components by estimated energy and ensuring
a consistent sign convention.

# Arguments
- `V::AbstractMatrix{<:Real}`  
  Demixing (unmixing) matrix of size `mxm` (e.g., returned by a joint diagonalizer).
- `W::AbstractMatrix{<:Real}`  
  Mixing matrix of size `mxm` (inverse of the true source mixing).

# Returns
- `B::Matrix{Float64}`  
  Separation matrix of size `mxm` such that `B * X` yields estimated source signals.
  Rows of `B` are ordered by descending source energy and each row is scaled to
  have a positive first coefficient (up to a small offset).

"""
function separate_sources(V::AbstractMatrix, W::AbstractMatrix)
    B = V' * W

    iW = inv(W)
    A = iW * V
    energies = vec(sum(abs2, A; dims=1))
    idx = sortperm(energies)

    B = B[idx[end:-1:1], :]

    first_row = B[1, :]
    signs = sign.(first_row .+ 0.1)
    B = Diagonal(signs) * B

   return B
end

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
    n,T = size(X)

    CM = cumulant_matrices(Matrix(X), n)
    V, _ = joint_diagonalization(CM, T, 2)
    
    #B = separate_sources(V, W)
    B = V' * W
    S = B * X
    iS = S'
    return hcat(time, iS)

end