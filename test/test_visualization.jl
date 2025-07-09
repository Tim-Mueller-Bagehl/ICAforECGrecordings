@testset "Visualization" begin
    # Write your tests here.
    @testset "plot_dataset" begin
        data = [0.0 1.0 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5] 
        plt = ICAforECGrecordings.plot_dataset(data) 
        @test !isnothing(plt)
    end
end