using LinearAlgebra
using Infinity
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

    N, T = size(X) 

    W = I(N)
    Y = X

    s_list = zeros(0)
    y_list = zeros(0)
    r_list = zeros(0)
    current_loss = ∞
    sign_change = false

    for i in maxiter
        psiY = score(Y)
        psidY_mean = score_der(psiY)

        g = gradient(Y, psiY)

        K = psidY_mean - diagm(g)
        signs = sign.(K)
        if n > 1
            sign_change = signs ≠ old_signs
        end
        old_signs = signs

        psidY_mean = psidY_mean .* signs

        G .= (g .- transpose(g)) ./ 2.0

        gradient_norm = max(max(abs.(G)))
        if gradient_norm < tol
            break
        end

        if n > 1
            append!(s_list, direction)
            y = G .- G_old
            append!(y_list, y)
            append!(r_list, 1. / sum(sum(direction .* y)))
            if length(s_list) > m
                deleteat!(s_list,1)
                deleteat!(y_list,1)
                deleteat!(r_list,1)
            end
        end
        G_old = G
        if sign_change
            s_list = zeros(0)
            y_list = zeros(0)
            r_list = zeros(0)
            current_loss = ∞
        end
        h = proj_hessian_approx(Y, psidY_mean, g)
        h = regularize_hessian(h, lambda_min)
        converged, new_Y, new_loss, alpha = line_search(Y, signs, direction, current_loss, ls_tries)

        if !converged
            direction = -G
            s_list = zeros(0)
            y_list = zeros(0)
            r_list = zeros(0)
            tmp, new_Y, new_loss, alpha = line_search(Y, signs, direction, current_loss, ls_tries)
        end
        direction = alpha * direction
        Y = new_Y
        W = expm(direction) .* W
        current_loss = new_loss
        if verbose
            @info iteration n, gradient norm = gradient_norm
        end
    end

end