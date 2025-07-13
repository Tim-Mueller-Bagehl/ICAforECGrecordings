include("StandardPicard.jl")
include("Picardo.jl")
"""
    function picard(X :: AbstractMatrix, params ::Dict)

Runs the Picard algorithm for ICA
#Arguments
- `X::AbstractMatrix` : Matrix that contains the signal that have to be unmixed
- `params::Dict` : Dictionary of optional Parameters
#Optional Arguments
- `m::Int` : Size for L-BFGS memory Typical values range between 3-15 default 7
- `maxiter::Int` : Maximal number of iterations for the algorithm default 100
- `mode::String` : choose between orthogonal picard ->"ortho"(deafult) and standard picard -> "standard"
- `tol::Float` : choose tolerance for stopping criterion default 1e-8
- `lambda_min::Float` : Constant used to regularize the Hessian approximation default 1e-2
- `ls_tries::Int` : Number of tries allowed for the backtracking line search default 10
- `whiten::Bool` : If true, the signals are whitened befor running ICA default true
- `verbose::Bool` : If true, prints the information about the algorithm default false
- `w_init::AbstractMatrix` : Initial rotation matrix for the algorithm default empty
- `python_defaults::Bool` : If true uses Python-compatible defaults for the algorithm default false
- `distribution::String` : Distribution used for the distribution-based ICA Possible values: "logistic"(default) or "logcosh"
- `renormalization::String` : Renormalization method used for the distribution-based ICA Possible values: "original" (default) or "pythonlike"
"""
function picard(X :: AbstractMatrix, params = Dict() ::Dict)

    if isempty(X)
        error("no signals Provided")
    end

    if length(size(X)) > 2
        error("Input signal must be two dimensional")
    end
    N,T = size(X)

    if N > T
        error("There are more signals than samples")
    end


    lowercase_parms = Dict(lowercase(key) => v for (key,v) in params)
    

    m = get(lowercase_parms,"m",7)
    maxiter = get(lowercase_parms,"maxiter",100)
    mode = get(lowercase_parms,"mode","ortho")
    tol = get(lowercase_parms,"tol",1e-8)
    lambda_min = get(lowercase_parms,"lambda_min",0.01)
    ls_tries = get(lowercase_parms,"ls_tries",10)
    whiten = get(lowercase_parms,"whiten",true)
    verbose = get(lowercase_parms,"verbose",false)
    centering = get(lowercase_parms,"centering",false)
    n_components = get(lowercase_parms,"pca",size(X,1))#whitening_mode
    w_init = get(lowercase_parms,"w_init",zeros(0))
    python_defaults = get(lowercase_parms,"python_defaults",false)
    distribution = get(lowercase_parms,"distribution","logistic")
    renormalization = get(lowercase_parms,"renormalization","original")
    if haskey(lowercase_parms,"pca")
        whitening_mode = "pca"
    else
        whitening_mode = "sph"
    end

    if python_defaults
        maxiter = 512
        tol = 1e-7
        m = 10
        centering = true
        distribution = "logcosh"
        renormalization = "pythonlike"
        whitening_mode = "pca"
        if verbose
            print("Using Python-compatible defaults")
        end
    end

    if isempty(w_init)
        if python_defaults
            if verbose
                print("Using random w_init")
            end
            w_init = random_init(n_components)
            w_init = sym_decorrelation(w_init)
        else
            w_init = I(n_components)
        end
    end

    if whiten == false && n_components != size(X,1)
        error("PCA works only if whiten = true")
    end

    if n_components != rank(X)
        @warn "Input matrix is of deficient rank. Please consider to reduce dimensionality (pca) prior to ICA"
    end

    if centering
        X_mean = mean(X,2)
        X = X .- X_mean
    end

    if whiten
        X_whitened , W_whitened = whitening_picard(X,whitening_mode,n_components)
    else
        X_whitened = X
        W_whitened = I(n_components)
    end
    if isempty(w_init)
        if python_defaults
            if verbose
                print("Using random w_init")
            end
            w_init = random_init(n_components)
        else
            w_init = I(n_components)
        end
    end
    X_whitened = w_init * X_whitened

    if mode == "ortho"
        Y , W_algo = picardo(X_whitened,m,maxiter,tol,lambda_min,ls_tries,verbose)
    elseif mode == "standard"
        Y , W_algo = picard_standard(X_whitened,m,maxiter,2,tol,lambda_min,ls_tries,verbose,distribution,renormalization)
    end

    W = W_algo * w_init * W_whitened
    return Y , W
end


function random_init(n_components :: Int)
    Q,R = qr(randn(n_components,n_components))
    return Matrix(Q)
end

function whitening_picard(Y,mode,n_components)

    R = (Y*transpose(Y)) / size(Y,2)
    U,D,_ = svd(R)
    #D = diagm(D)
    if mode == "pca"
        W = diagm(1 ./ sqrt.(D))*transpose(U)
        W = [start:n_components,:]
        for i = size(W,1)
            _,j = findmax(abs(W[i,:]))
            if W[i,j]<0
                W[i,:] = -W[i,:]
            end
        end
        Z = W * Y
    elseif mode == "sph"
        W = U* diagm(1 ./ sqrt.(D)) * transpose(U)
        Z = W * Y
    else
        error("$mode is not a valid mode. The options are pca and sph")
    end
    return Z,W
end

function laplace_rnd(shape ,mu :: Z= 0, b :: T= (1/sqrt(2))) where {T,Z<:Real}
    u = rand(shape) - 0.5
    y = mu - b * sign.(u) .* log(1 -2 *abs(u))
    return y  
end