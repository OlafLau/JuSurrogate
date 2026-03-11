using PhaseDiagram: DIS

# We test the polymer blend system AB3/AB3
# χN=19, f=0.4, α=0.14, C=1.0

@testset "free_energy.jl: DIS" begin
    # ϕ = 0.4
    _, config = init_model(DISPhase, "DIS.yml")
    # Evaluate at ϕ=0.4
    F, μ = free_energy(DISPhase, config)
    @test F ≈ 3.6134535528794847 atol=1e-6
    @test μ ≈ 0.6744577118259163 atol=1e-6
end

@testset "free_energy.jl: LAM" begin
    # ϕ = 0.4
    model, config = init_model(LAMPhase, "LAM.yml")
    # Evaluate at ϕ=0.4
    F, μ = free_energy(LAMPhase, config, model)
    @test F ≈ 3.3846565173210426 atol=1e-6 # ϕ=0.4
    @test μ ≈ 0.6940957577373614 atol=1e-4 # ϕ=0.4
    # Evaluate at ϕ=0.2
    # F, μ = free_energy(0.2, LAMPhase, config, model)
    # @test F ≈ 2.640834192563166 atol=1e-6 # ϕ=0.2
    # @test μ ≈ 0.2725506690831303 atol=1e-4 # ϕ=0.2
end

@testset "free_energy.jl: HEX" begin
    # Input is ϕ=0.4
    model, config = init_model(LAMPhase, "HEX.yml")
    # Evaluate at ϕ=0.2
    F, μ = free_energy(0.2, HEXPhase, config, model)
    @test F ≈ 2.6476634258093696 atol=1e-5 # ϕ=0.2
    @test μ ≈ 0.3280542449018094 atol=1e-4 # ϕ=0.2
end

@testset "free_energy.jl: DIS-DIS" begin
    ϕ₁ = 0.4
    ϕ₂ = 0.6
    ϕ₀ = 0.5
    m = MacrophaseModel("DIS.yml"; ϕ₀=ϕ₀)
    G, dG1, dG2 = gibbs_free_energy(ϕ₁, ϕ₂, m)
    χN = m.config1["Model"]["chiN"][1]
    fs = m.config1["Model"]["f"]
    f = fs[1]
    α = fs[1] * fs[end]
    tG, tdG1, tdG2 = DIS.G(ϕ₁, ϕ₂, ϕ₀, χN, f, α)
    @test G ≈ tG atol=1e-6
    @test dG1 ≈ tdG1 atol=1e-6
    @test dG2 ≈ tdG2 atol=1e-6
end

@testset "free_energy.jl: DIS-LAM" begin
    # DIS: ϕ=0.2, LAM ϕ=0.4
    ϕ₁ = 0.2
    ϕ₂ = 0.4
    ϕ₀ = 0.3
    m = MacrophaseModel("DIS.yml"; configfile2="LAM.yml", phase2=LAMPhase, ϕ₀=ϕ₀)
    G, dG1, dG2 = gibbs_free_energy(ϕ₁, ϕ₂, m)
    @test G ≈ 3.1147525366322153 atol=1e-6
    @test dG1 ≈ -0.19210901949156667 atol=1e-4
    @test dG2 ≈ 1.1293935170464424 atol=1e-4
end

@testset "free_energy.jl: LAM-LAM" begin
    # LAM: ϕ=0.2, LAM ϕ=0.4
    ϕ₁ = 0.2
    ϕ₂ = 0.4
    ϕ₀ = 0.3
    m = MacrophaseModel("LAM.yml"; phase1=LAMPhase, phase2=LAMPhase, ϕ₀=ϕ₀)
    G, dG1, dG2 = gibbs_free_energy(ϕ₁, ϕ₂, m)
    @test G ≈ 3.012745354942104 atol=1e-6
    @test dG1 ≈ -0.8861605651692256 atol=1e-4
    @test dG2 ≈ 0.6193576085958856 atol=1e-4
end