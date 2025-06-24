using ICAforECGrecordings
using Test

@testset "ICAforECGrecordings.jl" begin

    include("test_parser.jl")

    include("test_visualization.jl")

    include("test_preprocessing.jl")
end

