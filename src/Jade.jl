"""
    estimate_cumulant_matrices(X::AbstractMatrix) -> Matrix{E}

Compute the set of fourth-order cumulant matrices for the signals stored as rows of `X`.

# Arguments
- `X::AbstractMatrix{<:Number}`  
  An nxT data matrix, where each of the n rows is a sequence of T observations.

# Returns
- `CM::Matrix{E}` where `E = eltype(X)`  
  An nx(n*nbcm) matrix formed by horizontally concatenating `nbcm = n*(n+1)÷2` individual nxn cumulant blocks.  
  Each block corresponds to the fourth-order cumulant of the i-th and j-th signals (with i ≥ j).

"""
function estimate_cumulant_matrices(X::AbstractMatrix)
    E = eltype(X)

    n, T = size(X)
    nbcm = div(n * (n + 1), 2)
    CM = zeros(E, n, n * nbcm)
    R = one(E)*I(n)
    scale = fill(one(E)/T, n)

    idx = 1
    for i in 1:n, j in 1:i
        xi = @view X[i, :] # 1 x T
        xj = @view X[j, :] # 1 x T
        M = zeros(n, n)

        # Mij = (1/T) * [xi[t] * xj[t]] * X
        M = scale .* (xi .* xj)'
        Q = (M .* X) * X'

        # Diagonal Cumuant
        if i == j
            Q .-= R .+ 2*(R[:, i]*R[:, i]')
        # Off-Diagonal Cumulant
        else
            Q .-= (R[:, i]*R[:, j]' .+ R[:, j]*R[:, i]')
            Q .*= convert(E, sqrt(2))
        end
        CM[:, idx:idx+n-1] = Q
        idx += n
    end
    return CM
end

"""
    joint_diagonalization(CM_in::AbstractMatrix, T::Int; max_sweeps::Int=10)

Compute an approximate joint diagonalization of a set of `K` symmetric `mxm` matrices
that are concatenated horizontally in `CM_in`, using a Jacobi-style algorithm.

# Arguments
- `CM_in::AbstractMatrix{<:Real}`  
  Input matrix of size `m x (m*K)`.  Each contiguous block of `m` columns  
  corresponds to one symmetric matrix to be jointly diagonalized.
- `T::Int`  
  A scaling parameter (often the number of samples) used to set the rotation threshold.
- `max_sweeps::Int=10`  
  Maximum number of full sweeps over all index pairs before terminating.

# Returns
- `V::Matrix{Float64}`  
  An `mxm` orthogonal matrix whose columns jointly diagonalize the blocks.
- `CM::Matrix{Float64}`  
  The transformed version of `CM_in`, such that each `mxm` block is as close to  
  diagonal as possible under the rotation `V`.

"""
function joint_diagonalization(CM_in::AbstractMatrix, T::Int, max_sweeps::Int=10)
    CM = deepcopy(CM_in)
    m = size(CM, 1)
    K = size(CM, 2) ÷ m

    # 1) Init by diagonalizing first block
    V = eigen(Symmetric(CM[:, 1:m])).vectors # eigendecomposition Q1
    @views for k in 0:K-1 # rotate all blocks by V
        u = k*m + 1  : k*m + m
        @views CM[:, u] .*= V
    end
    CM .= V' * CM

    # 2) Joint Diagonalization by Jacobi
    threshold = 1/sqrt(T)/100
    sweep = 0
    updates = 0

    idxs = Tuple{Int,Int,AbstractRange{Int},AbstractRange{Int}}[]
    for p in 1:m-1, q in p+1:m
        push!(idxs, (p, q, p:m:m*K, q:m:m*K))
    end

    @inbounds while sweep < max_sweeps
        sweep += 1
        any_change = false

        for (p, q, Ip, Iq) in idxs
            @views begin
                g1 = CM[p, Ip] .- CM[q, Iq]
                g2 = CM[p, Iq] .+ CM[q, Ip]
            end
            gg11 = dot(g1,g1)
            gg22 = dot(g2,g2)
            gg12 = dot(g1,g2)
            ton = gg11 - gg22
            toff = gg12 + gg12
            theta = 0.5 * atan(toff, ton + sqrt(ton^2 + toff^2))

            if abs(theta) > threshold
                any_change = true
                updates += 1
                c = cos(theta); s = sin(theta)
                G = @views [c -s; s c]

                # update V on columns p,q
                @views V[:, [p,q]] .= V[:, [p,q]] * G

                # update CM on rows p,q
                @views CM[[p,q], :] .= G' * CM[[p,q], :]

                # update CM on columns Ip, Iq
                @views begin
                    tmpIp = copy(CM[:, Ip])
                    tmpIq = copy(CM[:, Iq])
                    CM[:, Ip] .= c .* tmpIp .+ s .* tmpIq
                    CM[:, Iq] .= -s .* tmpIp .+ c .* tmpIq
                end
            end
        end
        isnothing(any_change) && break
    end
    return V, CM
end

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
- `data::AbstractMatrix{<:Real}`  
  An `Nx(m+1)` matrix where the first column is a time vector of length `N`, and the
  remaining `m` columns are the observed mixed signals.

# Returns
- `result::Matrix{Float64}`  
  An `Nx(m+1)` matrix whose first column is the original time vector, and whose next
  `m` columns are the estimated source signals, ordered as they emerge from the
  JADE separation process.

"""
function jade(data::AbstractMatrix)
    data_w, W = whiten(data)

    time = data_w[:, 1]
    Xs = data_w[:, 2:end]
    X = Xs'
    _,T = size(X)
    CM = estimate_cumulant_matrices(Matrix(X))
    V, _ = joint_diagonalization(CM, T)
    
    #B = separate_sources(V, W)
    B = V' * W
    S = B * X
    iS = S'
    return hcat(time, iS)

end