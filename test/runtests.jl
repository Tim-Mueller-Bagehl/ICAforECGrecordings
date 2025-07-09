using ICAforECGrecordings
using Test
using Random: seed!
using Statistics

@testset "ICAforECGrecordings.jl" begin
    seed!(1234) 
    include("test_parser.jl")
    include("test_visualization.jl")
    include("test_preprocessing.jl")
    @testset "load_example_data" begin
        data = load_example_data()
        @test data isa Matrix{Float64}
    end

    include("test_shibbs.jl")
    


end

