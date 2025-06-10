using ICAforECGrecordings
using Test

@testset "parser" begin
    datpath = joinpath(@__DIR__, "..", "data", "foetal_ecg.dat")
    result = ReadDatasetFromDatFile(datpath)
    @test result isa Matrix{Float64}
end
