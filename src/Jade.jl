using LinearAlgebra

function estimate_cumulant_matrices(X::Matrix{Float64})
    n, T = size(X)
    nbcm = div(n * (n + 1), 2)
    CM = zeros(n, n * nbcm)
    R = I(n)
    scale = fill(1/T, n)

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
            Q .*= sqrt(2)
        end
        CM[:, idx:idx+n-1] = Q
        idx += n
    end
    return CM
end

function joint_diagonalization(CM_in::Matrix{Float64}, T::Int, max_sweeps::Int=10)
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

function separate_sources(V::Matrix{Float64}, W::Matrix{Float64})
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

function jade(data::Matrix{Float64})
    data_w, W = whiten(data)

    time = data_w[:, 1]
    Xs = data_w[:, 2:end]
    X = Xs'
    n,T = size(X)

    CM = estimate_cumulant_matrices(Matrix(X))
    V, _ = joint_diagonalization(CM, T)
    
    #B = separate_sources(V, W)
    B = V' * W
    S = B * X
    iS = S'
    return hcat(time, iS)
end
