using Test
using CubicHermiteSpline
using Polymer
using Scattering
using Polyorder: Polyorder, BB, SD, Config
using PhaseDiagram
using Random: Xoshiro

function model0()
    system = Polymer.AB_A_system()
    uc = UnitCell(1.0)
    lattice = BravaisLattice(uc)
    w = AuxiliaryField(zeros(1), lattice)
    nccscft = NoncyclicChainSCFT(system, w, 0.2; mde=OSF, init=:randn, rng=Xoshiro(3083))
    return PolyorderModel(nccscft)
end

function model1()
    system1 = Polymer.AB_A_system(; ϕAB=0.9)
    uc1 = UnitCell(1.0)
    lattice1 = BravaisLattice(uc1)
    w1 = AuxiliaryField(zeros(1), lattice1)
    nccscft1 = NoncyclicChainSCFT(system1, w1, 0.1; mde=OSF, init=:randn, rng=Xoshiro(3083))
    scftconfig = SCFTConfig(; max_iter=1000)
    config = Polyorder.Config(; scft=scftconfig)
    opt = VariableCellOpt(VariableCell(BB(1.0), SD(0.2)))
    return PolyorderModel(nccscft1, opt, config, false)
end

function model2()
    system2 = Polymer.AB_A_system(; ϕAB=0.9)
    uc2 = UnitCell(4.0)
    lattice2 = BravaisLattice(uc2)
    w2 = AuxiliaryField(zeros(64), lattice2)
    nccscft2 = NoncyclicChainSCFT(system2, w2, 0.01; mde=OSF, init=:randn, rng=Xoshiro(3083))
    scftconfig = SCFTConfig(; max_iter=1000)
    config = Config(; scft=scftconfig)
    opt = VariableCellOpt(VariableCell(BB(1.0), SD(0.2)))
    return PolyorderModel(nccscft2, opt, config, false)
end

@testset "phasemodel.jl: MacrophaseModel" begin
    mapm = MacrophaseModel(model0(), model0())
    @test mapm.phase1 == DISPhase
    @test mapm.model1.optimized == false
    @test mapm.model2.optimized == false
    @test isempty(mapm.kwargs) == true

    ϕc = ϕControlParameter(2)
    mpmodel = MacrophaseModel(model1(), model2(), ϕc; ϕ₁=0.2, ϕ₂=0.6, ϕ₀=0.4, cached=true)
    gibbs_free_energy(0.2, 0.6, mpmodel)
    gibbs_free_energy(0.1, 0.7, mpmodel)
    gibbs_free_energy(0.18, 0.62, mpmodel)
    m1 = PhaseDiagram.best_model(mpmodel, 1, 0.11)
    m2 = PhaseDiagram.best_model(mpmodel, 2, 0.88)
    @test getparam(m1, ϕc) == 0.1
    @test getparam(m2, ϕc) == 0.7

    @test length(mpmodel.model1_history) == 3
    PhaseDiagram.cache_model!(mpmodel, 1, m1)
    @test length(mpmodel.model1_history) == 3
    setparam!(m1, 0.11, ϕc)
    PhaseDiagram.cache_model!(mpmodel, 1, m1)
    @test length(mpmodel.model1_history) == 4
    @test [getparam(m, ϕc) for m in mpmodel.model1_history] == [0.2, 0.1, 0.18, 0.11]
end

@testset "phasemodel.jl: MacrophaseModel: LAMPhase" begin
    ϕc = ϕControlParameter(2)  # concentration of the second component, i.e. the homopolymer A.
    mapm = MacrophaseModel(model1(), model2(), ϕc; phase1=DISPhase, phase2=LAMPhase, cached=true)
    F1, μs1 = free_energy(0.1, 2, mapm)
    F2, μs2 = free_energy(0.2, 2, mapm)
    F3, μs3 = free_energy(0.3, 2, mapm)
    @test length(mapm.model2_history) == 3
    F4, μs4 = free_energy(0.2, 2, mapm)
    # Check for same parameter, no cache done!
    @test length(mapm.model2_history) == 3
    @test F4 ≈ F2 atol=1e-10
end

@testset "phasemodel.jl: SurrogateMacrophaseModel" begin
    ϕs, Fs, μs = rand(10), rand(10), rand(10)
    surrogate = CHSplineSurrogate(ϕs, Fs, μs)
    smapm = SurrogateMacrophaseModel(surrogate, surrogate)
    @test smapm.phase1 == DISPhase
    @test smapm.surrogate1.spl isa CubicHermiteSplineInterpolation
    @test smapm.surrogate2.spl isa CubicHermiteSplineInterpolation
    @test isempty(smapm.kwargs) == true
end

@testset "phasemodel.jl: MicrophaseModel" begin
    fc = fControlParameter(1, 1, fParam)  # f of the first block (A block) of the first component (AB diblock copolymer).
    mipm = MicrophaseModel(model1(), model2(), fc)
    @test mipm.ctrlparam == fc
    @test mipm.phase1 == DISPhase
    @test mipm.phase2 == LAMPhase
    @test mipm.model1.optimized == false
    @test mipm.model2.optimized == false
    @test isempty(mipm.kwargs) == true

    mipm = MicrophaseModel(model1(), model2(), fc; phase1=DISPhase, phase2=LAMPhase, cached=true)
    @test isempty(mipm.model1_history)
    @test isempty(mipm.model2_history)
    free_energy(0.1, 1, mipm)
    free_energy(0.2, 1, mipm)
    free_energy(0.3, 1, mipm)
    @test length(mipm.model1_history) == 3
    m1 = PhaseDiagram.best_model(mipm, 1, 0.18)
    @test PhaseDiagram.getparam(m1, fc) == 0.2
    m2 = PhaseDiagram.best_model(mipm, 1, 0.3)
    @test PhaseDiagram.getparam(m2, fc) == 0.3
    PhaseDiagram.cache_model!(mipm, 1, m2)
    @test length(mipm.model1_history) == 3
    setparam!(m1, 0.18, fc)
    PhaseDiagram.cache_model!(mipm, 1, m1)
    @test length(mipm.model1_history) == 4
    @test [getparam(m, fc) for m in mipm.model1_history] == [0.1, 0.2, 0.3, 0.18]

    F1, μs1 = free_energy(0.5, 2, mipm)
    F2, μs2 = free_energy(0.45, 2, mipm)
    F3, μs3 = free_energy(0.55, 2, mipm)
    @test length(mipm.model2_history) == 3
    F4, μs4 = free_energy(0.45, 2, mipm)
    # Check for same parameter, no cache done!
    @test length(mipm.model2_history) == 3
    @test F4 ≈ F2 atol=1e-10
end

@testset "phasemodel.jl: SurrogateMicrophaseModel" begin
    fc = fControlParameter(1, 1, fParam)
    ϕs, Fs, μs = rand(10), rand(10), rand(10)
    surrogate = CHSplineSurrogate(ϕs, Fs, μs)
    smipm = SurrogateMicrophaseModel(surrogate, surrogate, fc)
    @test smipm.ctrlparam == fc
    @test smipm.phase1 == DISPhase
    @test smipm.phase2 == LAMPhase
    @test smipm.surrogate1.spl isa CubicHermiteSplineInterpolation
    @test smipm.surrogate2.spl isa CubicHermiteSplineInterpolation
    @test isempty(smipm.kwargs) == true
end

nothing