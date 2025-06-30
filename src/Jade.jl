using LinearAlgebra

"""
jade(data::Matrix{Float64}, m::Int = size(data,2)-1) -> Matrix{Float64}

Performs blind source separation on already-whitened data using the JADE algorithm.

# Arguments
- `data`: A Tx(n+1) matrix where:
  - `data[:,1]` is the time vector (length T).
  - `data[:,2:end]` are n observed, whitened signals.
- `m`: Number of independent components to extract (default = 2).

# Returns
- `M1::Matrix{Float64}`: Tx2 matrix where
  - column 1 is `time`
  - column 2 is the first separated source signal
- `M2::Matrix{Float64}`: Tx2 matrix where
  - column 1 is `time`
  - column 2 is the second separated source signal
"""
function jade(data::Matrix{Float64}, m::Int = 2)
    # Separate time and signal data
    time = data[:,1]
    Xw   = data[:,2:end]             # T×n (time samples x channels)
    
    # Tranpose so we work with nxT (channels x samples)
    X = tranpose(Xw)
    n, T = size(X)
    @assert 1 ≤ m ≤ n "m must be between 1 and n"

    # 1) Build cumulant matrices for 4th-order statistics
    # flatten the symmetrix xumulant tensor into a set of n*(n+1)/2 matrices
    numCumulants = n * (n + 1) ÷ 2
    cumulants = zeros(n, n * numCumulants)
    I_n = I(n)
    scale = 1.0 / T #normalization constant

    idx = 1
    for i in 1:n
        xi = view(X, i, :)
        # Diagonal cumulant Q_{ii}
        M = zeros(n, n)
        for t in 1:T
            M .+= (xi[t]^2) * (X[:, t] * X[:, t]')
        end
        M .*= scale
        M .-= I_n .+ 2 * (I_n[:, i] * I_n[i, :])
        cumulants[:, idx:idx+n-1] = M
        idx += n

        # Off-diagonal cumulant Q_{ij} for i>j
        for j in 1:i-1
            xj = view(X, j, :)
            M = zeros(n, n)
            for t in 1:T
                M .+= (xi[t] * xj[t]) * (X[:, t] * X[:, t]')
            end
            M .*= scale * sqrt(2)
            M .-= (I_n[:, i] * I_n[j, :] + I_n[:, j] * I_n[i, :])
            cumulants[:, idx:idx+n-1] = M
            idx += n
        end
    end

    # 2) Joint diagonalization via Jacobi rotations
    V = I(n)
    thresh = 1/sqrt(T)/100
    encore = true

    while encore
        encore = false
        for p in 1:n-1, q in p+1:n
            Ip = collect(p:n:n*numCumulants)
            Iq = collect(q:n:n*numCumulants)
            g1 = cumulants[p, Ip] .- cumulants[q, Iq]
            g2 = cumulants[p, Iq] .+ cumulants[q, Ip]
            gg = [ dot(g1,g1)  dot(g1,g2);
                   dot(g2,g1)  dot(g2,g2) ]
            ton = gg[1,1] - gg[2,2]
            toff = gg[1,2] + gg[2,1]
            θ = 0.5 * atan(toff, ton + sqrt(ton^2 + toff^2))
            if abs(θ) > thresh
                encore = true
                c, s = cos(θ), sin(θ)
                G = [c -s; s c]
                V[:, [p,q]] .= V[:, [p,q]] * G
                cumulants[[p,q], :] .= G' * cumulants[[p,q], :]
                tmpIp, tmpIq = cumulants[:, Ip], cumulants[:, Iq]
                cumulants[:, Ip] .=  c .* tmpIp .+ s .* tmpIq
                cumulants[:, Iq] .= -s .* tmpIp .+ c .* tmpIq
            end
        end
    end

    # 3) Extract separated sources
    B = V'                     # n×n separating matrix
    S = B * X                  # n×T separated signals
    S = transpose(S)  # T×n back to time-major format

    # Create output matrices with time prefixed
    M1 = hcat(time, S[:, 1])
    M2 = hcat(time, S[:, 2])

    return M1, M2
end
