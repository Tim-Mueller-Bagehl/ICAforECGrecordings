"""
    cumulant_matrices(X::AbstractMatrix)
Blind separation of real signals with SHIBBS.
# Arguments
- `X::Matrix{Float64}`: Matrix that contains the signals that have to be unmixed, first Column is ignored. 
Returns 
CM=cumulant_matrices(X) a NxN*nbcm cumulant matrix.

"""
function cumulant_matrices(X::AbstractMatrix)
    N, T = size(X)
    CM = zeros(N, N*N)
    Rk = zeros(N, N)
    R = I(N)
    Xk = zeros(N, T)  
    xk = zeros(1, T) 
    
    for k in 1:N
        xk = X[k, :]
        Xk = X .* xk'
        Rk = Xk * Xk' / T - R
        Rk[k, k] -= 2
        CM[:, ((k - 1) * N + 1):(k * N)] .= Rk
    end

    return CM
end
"""
    joint_diagonalization(CM::AbstractMatrix, threshhold::Real)
Blind separation of real signals with SHIBBS.
# Arguments
- `CM::AbstractMatrix`: Matrix that contains the signals that have to be unmixed. 
-`threshhold::Real`: A threshhold for the for the rotation 
Returns 
V = joint_diagonalization(CM, threshhold) a NxN diagonalized matrix.

"""
function joint_diagonalization(CM::AbstractMatrix, threshhold::Real)
    N, _ = size(CM)
    V = Matrix{Float64}(I, N, N)
    nbrs = 1   
    sweep = 0    
    updates = 0 
    g	= zeros(2, N)
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
        for p = 1:N-1
            for q = (p+1):N
                Ip = p:N:N*N 
	            Iq = q:N:N*N
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
"""
    shibbs(X::Matrix{Float64},m::Int)
Blind separation of real signals with SHIBBS.
# Arguments
- `X::Matrix{Float64}`: whitend Matrix that contains the signals that have to be unmixed. 
Returns 
B=shibbsR(X) a nxn separating matrix such that S=B*X is an nxT
matrix of estimated source signals.

"""
function shibbs(X::Matrix{Float64})
    signals = X[:, 2:end]
    X = Matrix(signals')

    N, T = size(X)         
    threshhold = 0.01 / sqrt(T)
    MeanX   = mean(X, dims=2)
    X .-= MeanX 
    B = I(N)
    
    OneMoreStep = true
    outersweep = 0

    while OneMoreStep && outersweep < 1000

        outersweep += 1
        if(outersweep % 10 == 0)
            percent = outersweep / 10
            @info "outer loop done: " percent
        end

        CM = cumulant_matrices(X) # calculate cumulant matrices

        V = joint_diagonalization(CM, threshhold) # calculate diagonalization matrices

        RotSize = norm(V - I(N), 2) # calculate rotation size

        X = V' * X # Update X
        B = V' * B # Update B

        if RotSize < (N * threshhold) # check convergence
            OneMoreStep = false
        end  
    end

    # sort B By energy 
    A = pinv(B)
    vars = sum(A .* A, dims=2)
    keys = sortperm(vec(vars))
    B = B[keys, :]
    # fix signs of B to be non negativ 
    b = B[:, 1]
    signs = sign.(b)
    signs = [x == 0 ? 1 : x for x in signs]
    B = Diagonal(signs) * B

    return B
end