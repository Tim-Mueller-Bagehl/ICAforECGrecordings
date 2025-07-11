using ICAforECGrecordings: plot_dataset
using Test: @testset, @test, @test_throws

@testset "Visualization" begin
    # Write your tests here.
    @testset "plot_dataset" begin
        data = [0.0 1.0 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5] 
        plt = ICAforECGrecordings.plot_dataset(data) 
        @test !isnothing(plt)
    end

    @test "plot_dataset handles empty data" begin
        empty_data = Float64[]
        plt = ICAforECGrecordings.plot_dataset(empty_data)
        @test isnothing(plt)
    end
    @test "plot_dataset handles single signal" begin
        single_signal_data = [0.0 1.0; 1.0 2.0; 2.0 3.0]
        plt = ICAforECGrecordings.plot_dataset(single_signal_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles multiple signals" begin
        multi_signal_data = [0.0 1.0 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5]
        plt = ICAforECGrecordings.plot_dataset(multi_signal_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles non-numeric data" begin
        non_numeric_data = [0.0 "a" 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5]
        @test_throws MethodError ICAforECGrecordings.plot_dataset(non_numeric_data)
    end
    @test "plot_dataset handles NaN values" begin
        nan_data = [0.0 NaN 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5]
        plt = ICAforECGrecordings.plot_dataset(nan_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles Inf values" begin
        inf_data = [0.0 Inf 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5]
        plt = ICAforECGrecordings.plot_dataset(inf_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles large datasets" begin
        large_data = hcat(collect(0:9999), rand(10000, 2))
        plt = ICAforECGrecordings.plot_dataset(large_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles small datasets" begin
        small_data = [0.0 1.0; 1.0 2.0]
        plt = ICAforECGrecordings.plot_dataset(small_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles datasets with varying signal lengths" begin
        varying_length_data = [0.0 1.0; 1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0]
        plt = ICAforECGrecordings.plot_dataset(varying_length_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles datasets with zero signals" begin
        zero_signal_data = [0.0]
        plt = ICAforECGrecordings.plot_dataset(zero_signal_data)
        @test !isnothing(plt)
    end
    @test "plot_dataset handles sin cos functions" begin
        t = 0:0.01:10
        signal1 = sin.(2π .* t)
        signal2 = cos.(2π .* t)
        signal3 = sin.(2π .* t) .+ cos.(4π .* t)

        sin_data = hcat(t, signal1, signal2, signal3)        
        plt = ICAforECGrecordings.plot_dataset(sin_data)
        @test !isnothing(plt)
    end
end