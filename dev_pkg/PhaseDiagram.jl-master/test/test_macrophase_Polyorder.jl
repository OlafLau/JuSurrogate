@testset "macrophase.jl: MacrophaseModel" begin
    ϕ₀ = 0.7
    model = MacrophaseModel("DIS.yml"; ϕ₁=0.5, ϕ₂=0.9, ϕ₀=ϕ₀)

    ϕ₁, ϕ₂ = macrophase(OptimOptimizer(), model; lb=[0.1, ϕ₀], ub=[ϕ₀, 0.99])
    @test ϕ₁ ≈ 0.631055 atol=1e-3
    @test ϕ₂ ≈ 0.812722 atol=1e-3

    ϕ₁, ϕ₂ = macrophase(JSOOptimizer(), model; lb=[0.1, ϕ₀], ub=[ϕ₀, 0.99])
    @test ϕ₁ ≈ 0.631055 atol=1e-3
    @test ϕ₂ ≈ 0.812722 atol=1e-3
end

@testset "macrophase.jl: SurrogateMacrophaseModel" begin
    config_hex = load_config("HEX.yml")
    config_lam = load_config("LAM.yml")
    ϕs_hex, Fs_hex, μs_hex = readfile("free_energy_HEX.txt")
    ϕs_lam, Fs_lam, μs_lam = readfile("free_energy_LAM.txt")
    surrogate_hex = CHSplineSurrogate(ϕs_hex, Fs_hex, μs_hex, config_hex)
    surrogate_lam = CHSplineSurrogate(ϕs_lam, Fs_lam, μs_lam, config_lam)
    ϕ₀ = 0.183
    model = SurrogateMacrophaseModel(config_hex, surrogate_hex, surrogate_lam;
                                     config2=config_lam,
                                     ϕ₁=0.12, ϕ₂=0.24, ϕ₀=ϕ₀,
                                     phase1=HEXPhase, phase2=LAMPhase)

    ϕ₁, ϕ₂ = macrophase(OptimOptimizer(), model; lb=[0.02, ϕ₀], ub=[ϕ₀, 0.6])
    @test ϕ₁ ≈ 0.1758 atol=1e-3
    @test ϕ₂ ≈ 0.1904 atol=1e-3

    ϕ₁, ϕ₂ = macrophase(JSOOptimizer(), model; lb=[0.02, ϕ₀], ub=[ϕ₀, 0.6])
    @test ϕ₁ ≈ 0.1758 atol=1e-3
    @test ϕ₂ ≈ 0.1904 atol=1e-3

    ϕ₁, ϕ₂ = macrophase(RootsOptimizer(), model)
    @test ϕ₁ ≈ 0.1758 atol=1e-3
    @test ϕ₂ ≈ 0.1904 atol=1e-3
end