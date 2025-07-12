using ICAforECGrecordings
using Test
using Statistics: mean

@testset "whiten function" begin
    time = collect(1:10)
    X = [1.0 2.0;
         2.0 3.0;
         3.0 4.0;
         4.0 5.0;
         5.0 6.0;
         6.0 7.0;
         7.0 8.0;
         8.0 9.0;
         9.0 10.0;
         10.0 11.0]
    data = hcat(time, X)

    data_w, W = whiten(data)
    println("data_w type: ", typeof(data_w)) 
    
    time_white = data_w[:,1]
    X_white = data_w[:, 2:end]

    # Test time column
    @test isapprox(time, time_white)
    # Test mean
    @test all(abs.(mean(X_white, dims=1)) .< 1e-10)
    # Test size of return value
    @test size(W, 1) == size(X,2)
    @test size(W, 2) == size(X,2)
end