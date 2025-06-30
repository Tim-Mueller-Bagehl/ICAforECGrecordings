using ICAforECGrecordings
using Test
using Random
using Statistics

@testset "jade" begin
    Random.seed!(42)
    
    # create two non-gaussian signal
    T = 500
    time = collect(1:T)
    s1 = randn(T).^3
    s2 = randn(T).^3
    Strue = hcat(s1, s2)

    # mix it randomly
    A = [1.0 0.5; 0.2 1.2]
    Xmixed = hcat(time, Strue * A')

    # whitening
    Xw = whiten(Xmixed)

    # performing Jade
    M1, M2 = jade(Xw)

    # extract separated signals
    sep1 = M1[:, 2]
    sep2 = M2[:, 2]

    # caculate correlation
    c11 = abs(cor(sep1, s1)); c12 = abs(cor(sep1, s2))
    c21 = abs(cor(sep2, s1)); x22 = abs(cor(sep2, s2))

    # evaluate
    @test max(c11, c12) > 0.9
    @test max(c21, c22) > 0.9

end