@testset "parser" begin
    datpath = joinpath(@__DIR__, "..", "data", "foetal_ecg.dat")
    result = read_dataset_from_dat(datpath)
    @test result isa Matrix{Float64}
end
