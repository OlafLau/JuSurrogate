@testset "parameters.jl" begin
    # Polyorder configuration file
    config = load_config("DIS.yml")
    @test getparam(config, ϕParam) == 0.4
    @test getparam(config, χNParam) == 19.0
    @test getparam(config, fParam) == 0.4
    @test getparam(config, αParam) ≈ 0.14 atol=1e-15

    setparam!(config, 0.2, ϕParam)
    @test getparam(config, ϕParam) == 0.2
    @test config["Model"]["phi"][1] == 0.8
    setparam!(config, 20.0, χNParam)
    @test getparam(config, χNParam) == 20.0
    setparam!(config, 0.2, αParam)
    @test getparam(config, αParam) ≈ 0.2 atol=1e-15
end