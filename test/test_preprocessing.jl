@testset "Preprocessing" begin

    @testset "whiten() Test" begin
        data = [0.0 1.0 2.0; 1.0 1.5 2.5; 2.0 2.0 3.0; 3.0 2.5 3.5] 
        time = data[:, 1]
        data_white = whiten(data)

        #test: check the size
        @test size(data_white) == size(data)

        #test: check the time unmodified
        @test data_white[:, 1] ≈ time

        #test: check whether whitened data are centered for each row
        x_white = data_white[:, 2:end]
    end

end