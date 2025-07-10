
"""
    picardo(X::AbstractMatrix, m::Int, maxiter::Int, tol::Z, lambda_min <: A, ls_tries:: Int, verbose::Bool) where{Z,A <:Real}
Runs the Picard algorithm.
# Arguments
- `X::AbstractMatrix`: Matrix that contains the signals that have to be unmixed. 
- `m::Int` :  Size of L-BFGS's memory.
- `maxiter::Int` : Maximal number of iterations for the algorithm
- `tol::<:Real` : tolerance for the stopping criterion. Iterations stop when the norm of the projected gradient gets smaller than tol.
- `lambda_min::<:Real` : Constant used to regularize the Hessian approximation. The eigenvalues of the approximation that are below lambda_min are shifted to lambda_min.
- `ls_tries::Int` : Numer of tries allowed for the backtracking line-search. When that number is exceeded, the direction is thrown away and the gradient is used instead
- `verbose::Bool` : If true, prints the informations about the algorithm.

"""
function picardo(X::AbstractMatrix, m::Int, maxiter::Int, tol::Z, lambda_min :: A, ls_tries:: Int, verbose::Bool) where{Z,A <:Real}

    N, T = size(X) 

    W = I(N)
    Y = X

    s_list = zeros(0)
    y_list = zeros(0)
    r_list = zeros(0)
    current_loss = ∞
    sign_change = false

    for i  = 1:maxiter
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

function score(Y :: AbstractMatrix)
    return tanh.(Y)
end

function score_der(psiY :: AbstractMatrix)
    return -mean(psiY.^2,2) + 1
end
"""
    gradient(Y :: AbstractMatrix,psiY :: AbstractMatrix)
Compute the gradient of the current signals
"""
function gradient(Y :: AbstractMatrix,psiY :: AbstractMatrix)
    T = size(Y,2)
    return (psiY * transpose(T))/T
end
"""
    proj_hessian_approx(Y :: AbstractMatrix,psidY_mean :: AbstractMatrix,G :: AbstractMatrix)

Computes the projected Hessian approximation
"""
function proj_hessian_approx(Y :: AbstractMatrix,psidY_mean :: AbstractMatrix,G :: AbstractMatrix)
    N = size(Y,1)
    diagonal = psidY_mean * ones(1,N)
    off_diag = diag(G)
    off_diag = repeat(off_diag,1,N)
    return 0.5 * (diagonal + transpose(diagonal) - off_diag - transpose(off_diag))
end
"""
    regularize_hessian(h :: AbstractMatrix,l :: T) where {T<:Real}

Clips the eigenvalues of h to l
"""
function regularize_hessian(h :: AbstractMatrix,l :: T) where {T<:Real}
    mask = h .< l
    h[mask] .= l
    return h
end

"""
    solve_hessian(G,h)

Returns the inverse Hessian times G
"""
function solve_hessian(G,h)
    return G./h
end
"""
    loss(Y :: AbstractMatrix,signs :: AbstractMatrix)

Returns the loss function, evaluated for the current signals
"""
function loss(Y :: AbstractMatrix,signs :: AbstractMatrix)
    output = 0
    N,T = size(Y)
    for ii in eachrow(Y)
        y = Y[ii,:]
        s = signs[ii,:]
        output = output + s *(sum(abs(y) + log1p(exp(-2 * abs(y))))) / T
    end
    return output
end

function l_bfgs_direction(G,h,s_list :: AbstractVector,y_list :: AbstractVector,r_list :: AbstractVector)
    q = G
    a_list=Any[]
    for ii in eachindex(s_list)
        s = s_list[end - ii +1]
        y = y_list[end - ii +1]
        r = r_list[end - ii +1]
        alpha = r * sum(sum(s .* q))
        push!(a_list,alpha)
        q = q - alpha * y
    end
    z = solve_hessian(q,h)
    for ii in eachindex(s_list)
        s = s_list[ii]
        y = y_list[ii]
        r = r_list[ii]
        alpha = a_list[end - ii + 1]
        beta = r * sum(sum(y.*z))
        z = z + (alpha - beta)  * s
    end    
    return -z
end

"""
    line_search(Y :: AbstractMatrix,signs :: AbstractMatrix,direction :: AbstractMatrix,current_loss :: T,ls_tries :: Int) where {T<:Real}

Performs a backtracking line search, starting from Y and W, in the direction direction
"""
function line_search(Y :: AbstractMatrix,signs :: AbstractMatrix,direction :: AbstractMatrix,current_loss :: T,ls_tries :: Int) where {T<:Real}

    alpha = 1
    if current_loss == ∞
        current_loss = loss(Y,signs)
    end
    for ii = 1:ls_tries
        Y_new = expm(alpha * direction) * Y
        new_loss = loss(Y_new,signs)
        if new_loss < current_loss
            converged = true
            return converged , Y_new , new_loss, alpha
        end
        alpha = alpha /2
    end
    converged = false
    return converged , Y_new , new_loss, alpha
end
