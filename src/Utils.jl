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