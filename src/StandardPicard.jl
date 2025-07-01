using LinearAlgebra

function picard_standard(X, m, maxiter, precon, tol, lambda_min, ls_tries, verbose)
    N,T = size(X)
    W = I(N)
    Y = X
    s_list = zeros(0)
    y_list = zeros(0)
    r_list = zeros(0)
    current_loss = loss(Y,W)

    for n_top in 1:maxiter
        
        thY = tanh(Y/2.0)

        G = (thY * transpose(Y)) / T - I(N)

        g_norm = max(max(abs(G)))
        if g_norm < tol
            break
        end
        if n_top > 1
            append!(s_list,direction)#sketchy at best ein Fehler at worst direction noch nicht definiert
            y = G - G_old
            append!(y_list,y)
            append!(r_list,(1.0/sum(sum8direction .* y)))
            if length(s_list) > m
                s_list = s_list[2:end]
                y_list = y_list[2:end]
                r_list = r_list[2:end]
            end
        end
        G_old = G

        direction = l_bfgs_direction(Y,thY,G,s_list,y_list,r_list,precon,lambda_min)

        converged,new_Y,new_loss,direction = line_search(Y,W,direction,current_loss,ls_tries,verbose)

        if !converged
            direction = -G
            s_list = zeros(0)
            y_list = zeros(0)
            r_list = zeros(0)
            tmp,new_Y,new_loss,direction = line_search(Y,W,direction,current_loss,10,false)
        end
        Y= new_Y
        W = new_W
        current_loss = new_loss
        if verbose
            print("iteration $n_top graddient norm = $g_norm")
        end
    end
end


function loss(Y::Matrix{T},W::Matrix{T}) where T<:Real
    N=size(Y,1)
    loss=-log(det(W))

    for n in 1:N
        y=Y[n,:]
        loss = loss + mean(abs(y)+2. * log1p(exp(-abs(y))))
    end
    return loss
end

function line_search(Y,W,direction,current_loss,ls_tries,verbose)

    N = size(Y,1)
    projected_W = direction * W
    alpha=1
    for a in 1:ls_tries
        Y_new = (I(N) + alpha * direction) * Y
        W_new = W + alpha * projected_W
        new_loss = loss(Y_new,W_new)
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


function l_bfgs_direction(Y,thY,G,s_list,y_list,r_list,precon,lambda_min)
    q = G
    a_list=Any[]
    for ii in eachindex(s_list)
        s = s_list[end - ii +1]
        y = y_list[end - ii +1]
        r = r_list[end - ii +1]
        alpha = r*sum(sum(y .* q))
        push!(a_list,alpha)
        q = q- alpha * y
    end
    z = solve_hessian(q,Y,thY,precon,lambda_min)
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

function solve_hessian(G,Y,thY,precon,lambda_min)
    N , T =size(Y)

    psidY = (- thY.^2 + 1.0)/2.0

    Y_squared = Y.^2
    if precon == 2
        a = (psidY * transpose(Y_squared)) / T
    elseif precon == 1
        sigma2 = mean(Y_squared,2)
        psidY_mean = mean(psidY,2)
        a = psidY * sigma2
        diagonal_term = mean(mean(Y_squared .* psidY)) +1.0
        a[1:(N+1):N*N] .= diagonal_term
    else
        error("precon should be 1 or 2")
    end
    eigenvalues = 0.5 * (a + transpose(a) - sqrt((a-transpose(a)).^2 +4.0))

    problematic_locs = eigenvalues .< lambda_min
    problematic_locs[1:(N+1):N*N] .= false
    indices = findall(problematic_locs)
    i_pb = indices[1][1]
    j_pb=indices[1][2]
    a[i_pb,j_pb] += lambda_min - eigenvalues[i_pb,j_pb]

    return (G .* transpose(a) - transpose(G))./ (a.* transpose(a) - 1.0)

end


