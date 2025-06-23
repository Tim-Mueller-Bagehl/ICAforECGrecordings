"""
    picard(X::Matrix{Float64}, m::Int, maxiter::Int, tol::Float32, lambda_min::Float32, ls_tries, verbose::Bool)
Runs the Picard algorithm.
# Arguments
- `X::Matrix{Float64}`: Matrix that contains the signals that have to be unmixed. 
- `m::Int` :  Size of L-BFGS's memory.
- `maxiter::Int` : Maximal number of iterations for the algorithm
- `tol::Float32` : tolerance for the stopping criterion. Iterations stop when the norm of the projected gradient gets smaller than tol.
- `lambda_min::Float32` : Constant used to regularize the Hessian approximation. The eigenvalues of the approximation that are below lambda_min are shifted to lambda_min.
- `ls_tries` : Number of tries allowed for the backtracking line-search. When that number is exceeded, the direction is thrown away and the gradient is used instead
- `verbose::Bool` : If true, prints the informations about the algorithm.

"""
function picard(X::Matrix{Float64}, m::Int, maxiter::Int, tol::Float32, lambda_min::Float32, ls_tries, verbose::Bool)