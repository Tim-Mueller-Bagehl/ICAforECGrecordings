using ICAforECGrecordings: plot_dataset
using Test: @testset, @test, @test_throws

@testset "Visualization" begin
    # Write your tests here.
    @testset "plot_dataset" begin
        data = [0.0 1.0 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5] 
        plt = ICAforECGrecordings.plot_dataset(data) 
        @test !isnothing(plt)
    end

    @testset "plot_dataset with empty data" begin
        data = zeros(Float64, 0, 2) 
        plt = ICAforECGrecordings.plot_dataset(data)
        @test plt === nothing
    end

    @testset "plot_dataset with two signals" begin
        data = [0.0 1.0 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5]
        plt = ICAforECGrecordings.plot_dataset(data)
        @test !isnothing(plt)
    end
    @testset "plot_dataset with multiple signals" begin
        data = [0.0 1.0 2.0 3.0; 1.0 1.5 2.5 3.5; 2.0 2.0 3.0 4.0; 3.0 2.5 3.5 4.5]
        plt = ICAforECGrecordings.plot_dataset(data)
        @test !isnothing(plt)
    end



end