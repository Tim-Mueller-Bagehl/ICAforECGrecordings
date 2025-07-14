using ICAforECGrecordings
using Test
@testset "jade function" begin
    
    time = 0:pi/100:4*pi
    noise = rand(length(time)) .- 0.5   
    signal = sin.(time) 
    square_wave = squarewave1.(time)

    data = hcat(time, signal, noise)
    signals = solve(JadeSeperator(), data)
    @test size(signals) == size(data)
    @test isapprox(signals[:, 2], signal, atol=1e-5)
    @test isapprox(signals[:, 3], noise, atol=1e-5)

end
