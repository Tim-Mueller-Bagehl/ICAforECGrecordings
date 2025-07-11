"""
    shibbs(X::AbstractMatrix,m::Int)
Blind separation of real signals with SHIBBS.
# Arguments
- `X::AbstractMatrix`:  Matrix that contains the signals that have to be unmixed. 
- `m=nothing`: amount of signals that should be extracted from X. default same amount of singals as source
Returns 
B=shibbsR(X) a mxn separating matrix such that S=B*X is an mxT
matrix of estimated source signals.

"""
function shibbs(X::AbstractMatrix,m=nothing)

    signals = X[:, 2:end]
    X = Matrix(signals')
    N, T = size(X)   
    m = isnothing(m) ? N : m

    threshhold = 0.01 / sqrt(T)
    MeanX   = mean(X, dims=2)
    X .-= MeanX 
    C = (X * X') / T                    
    eig = eigen(C)                      

    D = eig.values                      
    U = eig.vectors                     

    k = sortperm(D)                         
    puiss = D[k]                        

    rangeW = N-m+1:N              
    scales = sqrt.(puiss[rangeW])       

    B = Diagonal(1.0 ./ scales) * U[:, k[rangeW]]'
    X = B * X  
    
    OneMoreStep = true
    outersweep = 0

    while OneMoreStep && outersweep < 1000

        outersweep += 1

        CM = cumulant_matrices(X,m) # calculate cumulant matrices

        V = joint_diagonalization(CM, threshhold,m) # calculate diagonalization matrices

        RotSize = norm(V - I(m), 2) # calculate rotation size

        X = V' * X # Update X
        B = V' * B # Update B

        if RotSize < (m * threshhold) # check convergence
            OneMoreStep = false
        end  
    end

    # sort B By energy 
    A = pinv(B)
    vars = sum(A .* A, dims=2)
    keys = sortperm(vec(vars))
    B = B[:, keys]
    # fix signs of B to be non negativ 
    b = B[:, 1]
    signs = sign.(b)
    signs = [x == 0 ? 1 : x for x in signs]
    B = Diagonal(signs) * B

    return B
end