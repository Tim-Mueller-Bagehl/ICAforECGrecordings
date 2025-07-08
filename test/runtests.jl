using ICAforECGrecordings
using Test

@testset "ICAforECGrecordings.jl" begin

    include("test_parser.jl")

    include("test_visualization.jl")

    include("test_preprocessing.jl")

    @testset "load_example_data" begin
        data = load_example_data()
        @test data isa Matrix{Float64}
    end
    


end

