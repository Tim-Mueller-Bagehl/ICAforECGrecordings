"""
    extract_sources(V::AbstractMatrix, W::AbstractMatrix) -> B::Matrix{Float64}

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
function extract_sources(V::AbstractMatrix, W::AbstractMatrix)
    B = V' * W

    iW = inv(W)
    A = iW * V
    energies = vec(sum(abs2, A; dims=1))
    idx = sortperm(energies)

    B = B[idx[end:-1:1], :]

    first_row = B[1, :]
    signs = sign.(first_row .+ 0.1)
    B = Diagonal(signs) * B

   return B[1:2, :]
end

"""
"""
function joint_diagonalization_new(CM_in::AbstractMatrix, T::Int, max_sweeps::Int=10)
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
    return V
end

"""
    fourth_order_cumulant_matrices(X::AbstractMatrix) -> Matrix{E}

Compute the set of fourth-order cumulant matrices for the signals stored as rows of `X`.

# Arguments
- `X::AbstractMatrix{<:Number}`  
  An nxT data matrix, where each of the n rows is a sequence of T observations.

# Returns
- `CM::Matrix{E}` where `E = eltype(X)`  
  An nx(n*nbcm) matrix formed by horizontally concatenating `nbcm = n*(n+1)÷2` individual nxn cumulant blocks.  
  Each block corresponds to the fourth-order cumulant of the i-th and j-th signals (with i ≥ j).

"""
function fourth_order_cumulant_matrices(X::AbstractMatrix)
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
    cumulant_matrices(X::AbstractMatrix, m::Number)
Blind separation of real signals with SHIBBS.
# Arguments
- `X::AbstractMatrix`: Matrix that contains the signals that have to be unmixed. 
- `m::Number`: The number of signals that should be extracted from X.
Returns 
CM=cumulant_matrices(X, m) a NxN*nbcm cumulant matrix.

"""
function cumulant_matrices(X::AbstractMatrix, m::Number)
    N, T = size(X)
    CM = zeros(m, m*m)
    Rk = zeros(m, m)
    R = I(m)
    Uns = ones(Int, m)  
    Xk = zeros(m, T)  
    xk = zeros(1, T) 
    
    for k in 1:m
        xk = X[k, :]
        Xk = X .* xk[Uns, :]  
        Rk = Xk * Xk' / T - R
        Rk[k, k] -= 2
        CM[:, ((k - 1) * m + 1):(k * N)] .= Rk
    end

    return CM
end

"""
    joint_diagonalization(CM::AbstractMatrix, threshhold::Real, m::Number)
Blind separation of real signals with SHIBBS.
# Arguments
- `CM::AbstractMatrix`: Matrix that contains the signals that have to be unmixed. 
-`threshhold::Real`: A threshhold for the for the rotation
- `m::Number`: The number of signals that should be extracted from CM. 
Returns 
V = joint_diagonalization(CM, threshhold, m) a NxN diagonalized matrix.

"""
function joint_diagonalization(CM::AbstractMatrix, threshhold::Real, m::Number)
    
    V = Matrix{Float64}(I, m, m)
    nbrs = 1   
    sweep = 0    
    updates = 0 
    g	= zeros(2, m)
    gg	= zeros(2, 2)
    G	= zeros(2, 2)
    c = 0.0
    s = 0.0
    ton	= 0.0
    toff = 0.0
    theta = 0.0
    repeat = true
    
    while repeat && sweep < 5000
        sweep += 1
        nbrs = 0
        for p = 1:m-1
            for q = (p+1):m
                Ip = p:m:m*m 
	            Iq = q:m:m*m
	            g = [CM[p,Ip] .- CM[q,Iq]; CM[p,Iq] .+ CM[q,Ip]]
	            gg	= g * g' 
	            ton = gg[1,1] - gg[2,2]
                toff = gg[1,2] + gg[2,1]
 	            theta	= 0.5 * atan(toff, ton + sqrt(ton^2 + toff^2))
	            if abs(theta) > threshhold	
                    nbrs += 1
	                c = cos(theta) 
	                s = sin(theta)
	                G = [c -s; s c]
	  
	                pair = [p, q] 
	                V[:, pair] = V[:, pair] * G
	                CM[pair, :] .= G' * CM[pair, :]
	                CM[:, [Ip; Iq]] = [c * CM[:, Ip] + s * CM[:, Iq] -s * CM[:, Ip] + c * CM[:, Iq]]
	            end
            end
        end
        updates += nbrs
        if nbrs == 0
            break
        end
    end
    return V
end