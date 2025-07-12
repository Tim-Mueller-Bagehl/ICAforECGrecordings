using ICAforECGrecordings
using Test
using LinearAlgebra: div, qr, Diagonal, I, norm, inv, adjoint, diag
using Random

@testset "extract_sources function" begin
    Random.seed!(123)

    m = 4
    n = 1000

    W = randn(m, m)
    V = Matrix(qr(randn(m, m)).Q)

    B = extract_sources(V,W)

    # Test size (2 x m)
    @test size(B) == (2, m)

    # Test row norm which should be distinctively bigger than 0
    norms = norm.(eachrow(B))
    @test all(norms .> 0.1)
end

@testset "joint_diagonalization_new" begin
    Random.seed!(123)

    m = 4
    K = 6
    T = 1000

    A = randn(m, m)
    # Use QR Decomposition for creating orthogonal matrix
    Q, _ = qr(A)
    Q = Matrix(Q)

    # Create multiple symmetric matrix for joint diagonalization
    CM_blocks = Matrix{Float64}[]
    for _ in 1:K
        D = Diagonal(randn(m))
        M = Q * D * Q'  # similarity transform M = QDQ' -> every M can be diagonalized by Q
        push!(CM_blocks, M)
    end

    CM = hcat(CM_blocks...)
    V = joint_diagonalization_new(CM, T)

    # Test the size of return value
    @test size(V) == (m, m)

    # Test orthonality
    @test isapprox(V' * V, I(m), atol=1e-10)

    #Test diagonality
    for k in 0:K-1  # check all matrices/blocks
        idx = k*m + 1 : (k+1)*m
        M_diag = V' * CM[:, idx] * V
        off_diag = M_diag .- Diagonal(diag(M_diag))
        rel_off = norm(off_diag) / norm(M_diag)
        @test rel_off < 0.05
    end

end

@testset "fourth_order_cumulant_matrices function" begin
    n, T = 3, 1000
    X = randn(n, T)
    CM = fourth_order_cumulant_matrices(X)
    nbcm = div(n * (n+1), 2)

    # Test the size of return value
    @test size(CM) == (n, n * nbcm)

    # Test the type or return value
    @test eltype(CM) == Float64

    # Test the symmetry of each cumulant matrix
    for k in 0:(nbcm-1)
        Q = CM[:, k*n+1:(k+1)*n]
        @test isapprox(Q, Q', atol=1e-10)
    end
end