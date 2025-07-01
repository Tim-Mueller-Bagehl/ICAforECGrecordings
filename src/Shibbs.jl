using LinearAlgebra
"""
    shibbs(X::Matrix{Float64},m::Int)
Blind separation of real signals with SHIBBS.
# Arguments
- `X::Matrix{Float64}`: Matrix that contains the signals that have to be unmixed. 
Returns 
B=shibbsR(X) a nxn separating matrix such that S=B*X is an nxT
matrix of estimated source signals.

"""
function shibbs(X::Matrix{Float64})
    signals = X[:, 2:end]
    X = Matrix(signals')
    N, T = size(X)         
    seuil = 0.01 / sqrt(T)
    MeanX   = mean(X, dims=2)
    X .-= MeanX 

    nbcm = N 
    CM = zeros(N, N*nbcm)
    Rk = zeros(N, N)
    R = I(N)
    Xk = zeros(N, T)  
    xk = zeros(1, T) 
    B = I(N)

    OneMoreStep = true
    outersweep = 0

    while OneMoreStep
        Range = 1:N
        for k = 1:N
            xk = X[k, :]
            xk = vec(X[k, :])
            Xk = X .* reshape(xk, 1, :) 
            Rk[k,k] -= 2
            CM[:, Range] = Rk   
            Range = Range .+ N
        end
    
        V = Matrix{Float64}(I, N, N)
        nbrs = 1   
        sweep = 0    
        updates = 0 
        g	= zeros(2, nbcm)
        gg	= zeros(2, 2)
        G	= zeros(2, 2)
        c = 0.0
        s = 0.0
        ton	= 0.0
        toff = 0.0
        theta = 0.0
        while true
            sweep += 1
            nbrs = 0
            for p = 1:N-1
                for q = (p+1):N
                    Ip = p:N:N*nbcm 
	                Iq = q:N:N*nbcm

	                g = [CM[p,Ip] .- CM[q,Iq]; CM[p,Iq] .+ CM[q,Ip]]
	                gg	= g * g' 
	                ton = gg[1,1] - gg[2,2]
                    toff = gg[1,2] + gg[2,1]
 	                theta	= 0.5 * atan(toff, ton + sqrt(ton^2 + toff^2))
	                if abs(theta) > seuil	
                        nbrs += 1
	                    c = cos(theta) 
	                    s = sin(theta)
	                    G = [c -s; s c]
	  
	                    pair = [p, q] 
	                    V[:, pair] = V[:, pair] * G
	                    CM[pair, :] .= G' * CM[pair, :]
	                    CM[:, vcat(Ip, Iq)] .= hcat(
                            c .* CM[:, Ip] .+ s .* CM[:, Iq],
                            -s .* CM[:, Ip] .+ c .* CM[:, Iq]
                        )
	                end
                end
            end

            updates += nbrs
            if nbrs == 0
                break
            end
            if sweep > 1000
                break
            end
        end
        outersweep += sweep
        RotSize = norm(V - I(N), 2)

        X = V' * X 
        B = V' * B
        if RotSize < (N * seuil)
            OneMoreStep = false
        end
        
        if outersweep > 1000 * 1000
                OneMoreStep = false
        end
    end

    A = pinv(B)
    vars = sum(A .* A, dims=2)
    keys = sortperm(vec(vars))
    B = B[keys, :]
    B = B[end:-1:1, :]

    b = B[:, 1]
    signs = sign.(b)
    signs = [x == 0 ? 1 : x for x in signs]
    B = Diagonal(signs) * B

    return B
end