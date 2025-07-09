using LinearAlgebra


"""
    picard_standard(X :: AbstractMatrix, m :: Int, maxiter :: Int, precon, tol :: Float32, lambda_min :: Float32, ls_tries :: Int, verbose :: Bool)

 Runs the Picard algorithm detailed in: https://arxiv.org/abs/1706.08171

# Arguments
- `X::AbstractMatrix`: Matrix that contains the signals that have to be unmixed. 
- `m::Int` :  Size of L-BFGS's memory.
- `maxiter::Int` : Maximal number of iterations for the algorithm
- `precon::Int` : Choose wich Hessian approximation is Used. 1 -> H1 2->H2 H2 is more costly but can greatly accelerate convergence
- `tol::Float32` : tolerance for the stopping criterion. Iterations stop when the norm of the projected gradient gets smaller than tol.
- `lambda_min::Float32` : Constant used to regularize the Hessian approximation. The eigenvalues of the approximation that are below lambda_min are shifted to lambda_min.
- `ls_tries` : Number of tries allowed for the backtracking line-search. When that number is exceeded, the direction is thrown away and the gradient is used instead
- `verbose::Bool` : If true, prints the informations about the algorithm.
- `distribution::String` : logistic or logcosh. The distribution to use for the score function
- `renormalization::String` : original or pythonlike. The method for Hessian renormalization
"""
function picard_standard(X :: AbstractMatrix, m :: Int, maxiter :: Int, precon :: Int, tol :: Float32, lambda_min :: Float32, ls_tries :: Int, verbose :: Bool, distribution ::String = "logistic",renormalization ::String = "original")
    @doc"""Init"""
    N,T = size(X)
    W = I(N)
    Y = X
    s_list = zeros(0)
    y_list = zeros(0)
    r_list = zeros(0)
    current_loss = loss(Y,W,distribution)

    for n_top = 1:maxiter
        @doc"""computes the score function"""
        if distribution == "logistic"
            thY = tanh(Y/2.0)
        else
            thY = tanh(Y)
        end
        @doc"""compute the relative gradient"""
        G = (thY * transpose(Y)) / T - I(N)
        @doc"""stopping criterion """

        g_norm = max(max(abs(G)))
        if g_norm < tol
            break
        end
        @doc"""update the memory """
        if n_top > 1
            append!(s_list,direction)
            y = G - G_old
            append!(y_list,y)
            append!(r_list,(1/sum(sum(direction .* y))))
            if length(s_list) > m
                s_list = s_list[2:end]
                y_list = y_list[2:end]
                r_list = r_list[2:end]
            end
        end
        G_old = G

        @doc"""Find the L-BFGS direction """
        direction = l_bfgs_direction(Y,thY,G,s_list,y_list,r_list,precon,lambda_min)
        @doc"""Do a line search in that direction """
        converged,new_Y,new_loss,direction = line_search(Y,W,direction,current_loss,ls_tries,verbose,distribution)

        if !converged
            direction = -G
            s_list = zeros(0)
            y_list = zeros(0)
            r_list = zeros(0)
            tmp,new_Y,new_W,new_loss,direction = line_search(Y,W,direction,current_loss,10,false,distribution)
        end
        Y= new_Y
        W = new_W
        current_loss = new_loss
        if verbose
            print("iteration $n_top graddient norm = $g_norm loss = $current_loss")
        end
    end
    return Y,W
end

"""
    loss(Y::AbstractMatrix,W::AbstractMatrix) 

compute the loss function for Y and W
"""
function loss(Y::AbstractMatrix,W::AbstractMatrix,distribution::String) 
    
    N=size(Y,1)
    loss=-log(det(W))
    if distribution == "logistic"
        for n = 1:N
            y=Y[n,:]
            loss = loss + mean(abs(y)+2 * log1p(exp(-abs(y))))
        end
        return loss
    else
        for k = 1:N
            y = Y[k,:]
            loss = loss + mean(abs(y) + 2 *log1p(exp(-abs(y))))
        end
        return loss
    end
end

"""
    line_search(Y :: AbstractMatrix,W :: AbstractMatrix,direction :: AbstractMatrix,current_loss :: Number,ls_tries :: Int,verbose :: Bool,distribution :: String)

Performs a backtracking line search, starting from Y and W, in the direction direction
"""
function line_search(Y :: AbstractMatrix,W :: AbstractMatrix,direction :: AbstractMatrix,current_loss :: Number,ls_tries :: Int,verbose :: Bool,distribution :: String)

    N = size(Y,1)
    projected_W = direction * W
    alpha=1
    for a = 1:ls_tries
        Y_new = (I(N) + alpha * direction) * Y
        W_new = W + alpha * projected_W
        new_loss = loss(Y_new,W_new,distribution)

        if new_loss < current_loss
            convergend = true
            rel_step = alpha * direction
            return convergend,Y_new,W_new,new_loss,rel_step
        end
        alpha = alpha/2
    end
    if verbose
        print("line search failed, falling back to gradient")
    end
    converged = false
    rel_step = alpha * direction
    return converged,Y_new,W_new,new_loss,rel_step
end


function l_bfgs_direction(Y :: AbstractMatrix,thY :: Number,G :: AbstractMatrix,s_list :: AbstractVector,y_list :: AbstractVector,r_list :: AbstractVector,precon :: Int,lambda_min :: Float32,distribution :: String,renormalization :: String)
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
    z = solve_hessian(q,Y,thY,precon,lambda_min,distribution,renormalization)
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

function regularize_hessian_pythonlike(a :: AbstractMatrix,lm :: Float32)
    N = size(a,1)

    e = 0.5*(a + transpose(a) - sqrt((a - transpose(a)).^2 + 4))

    mask = e .< lam
    mask[1:(N+1):N*N] .= false
    CartesianIndices =findall(mask)
    if !isempty(CartesianIndices)
        idx = LinearIndices(mask)[CartesianIndices]
        a[idx] = a[idx] + (lm - e[idx])
    end
    return a
    
end

function solve_hessian(G :: AbstractMatrix,Y :: AbstractMatrix,thY :: Number,precon :: Int,lambda_min::Float32,distribution :: String,renormalization :: String)
    N , T =size(Y)

    if distribution == "logistic"
        psidY = (- thY.^2 + 1)/2
    else
        psidY = 1- thY.^2
    end

    Y_squared = Y.^2
    if precon == 2
        a = (psidY * transpose(Y_squared)) / T
    elseif precon == 1
        sigma2 = mean(Y_squared,2)
        psidY_mean = mean(psidY,2)
        a = psidY_mean * transpose(sigma2)
        diagonal_term = mean(mean(Y_squared .* psidY)) +1
        a[1:(N+1):N*N] .= diagonal_term
    else
        error("precon should be 1 or 2")
    end

    if renormalization == "original"

        eigenvalues = 0.5 * (a + transpose(a) - sqrt((a-transpose(a)).^2 +4))

        problematic_locs = eigenvalues .< lambda_min
        problematic_locs[1:(N+1):N*N] .= false
        indices = findall(problematic_locs)
        i_pb = indices[1][1]
        j_pb=indices[1][2]
        a[i_pb,j_pb] = a[i_pb,j_pb] + lambda_min - eigenvalues[i_pb,j_pb]
    else
        a = regularize_hessian_pythonlike(a,lambda_min)
    end
    return (G .* transpose(a) - transpose(G))./ (a.* transpose(a) - 1)

end


