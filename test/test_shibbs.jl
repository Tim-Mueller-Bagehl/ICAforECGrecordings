# Function to compare two matrices up to sign and permutation
function is_equivalent_up_to_perm_sign(S_est, S_true; atol=1e-2)
    N, T = size(S_true)
    S_true_norm = S_true ./ std(S_true, dims=2)
    S_est_norm = S_est ./ std(S_est, dims=2)
    
    used = falses(N)
    for i in 1:N
        found = false
        for j in 1:N
            if !used[j]
                corr = abs(dot(S_true_norm[i, :], S_est_norm[j, :]) / T)
                if corr > 0.99
                    used[j] = true
                    found = true
                    break
                end
            end
        end
        if !found
            return false
        end
    end
    return true
end

@testset "shibbs" begin

    N = 3
    T = 20

    S = [sin.(2π .* (1:T) .* f ./ T) for f in [3.0, 7.0, 11.0]]
    S = reduce(vcat, S)
    S = reshape(S, N, T)

    A = randn(N, N)
    X = A * S

     X_dummy = hcat(zeros(N), X)

    B = shibbs(X_dummy)
    S_est = B * X


    @test is_equivalent_up_to_perm_sign(S_est, S), "Shibbs failed to recover the original sources up to permutation and sign."
    @test size(S_est) == size(S), "Shibbs output size does not match input size."
    @test all(ismissing, S_est) == false, "Shibbs output contains missing values."
    @test all(isnan, S_est) == false, "Shibbs output contains NaN values."
    @test all(isfinite, S_est) == true, "Shibbs output contains infinite values."   
end