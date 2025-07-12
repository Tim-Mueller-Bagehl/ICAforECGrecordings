using ICAforECGrecordings
using Test
using Random
using LinearAlgebra: adjoint, transpose, isapprox, randn
using Statistics: cor

#=
@testset "jade function" begin
    Random.seed!(123)

    m = 4
    n = 1000
    time = (collect(1:n))
    freq1, freq2 = 0.05, 0.07

    s1 = sin.(2*pi * freq1 * time)
    s2 = sin.(2*pi * freq2 * time)
    sources = [s1'; s2']

    W = randn(m, m)

    X = W * vcat(sources, randn(m-2, n))
    
    data_w = hcat(time, X')
    output = jade(data_w, W)

    # Test the size of return value
    @test size(output) == (n, 3)
    # Test the time column
    @test isapprox(output[:,1], time)

    # Test whether signals are Independent
    recovered = output[:, 2:3]
    norms = norm.(eachcol(recovered))
    @test all(norms .> 0.1)
    corr = cor(recovered[:, 1], recovered[:, 2])
    @test abs(corr) < 0.3
end
=#
