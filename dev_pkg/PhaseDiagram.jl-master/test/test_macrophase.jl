using Test
using Optim
using Roots
using Polymer
using Scattering
using Polyorder
using PhaseDiagram
import PhaseDiagram: DIS

function starAB3(; fA=0.4, fB1=0.2, fB2=0.2, fB3=0.2)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb0 = BranchPoint(:EB0)
    A = PolymerBlock(:A, sA, fA, FreeEnd(:A), eb0)
    B1 = PolymerBlock(:B1, sB, fB1, FreeEnd(:B1), eb0)
    B2 = PolymerBlock(:B2, sB, fB2, FreeEnd(:B2), eb0)
    B3 = PolymerBlock(:B3, sB, fB3, FreeEnd(:B3), eb0)
    return BlockCopolymer(:AB3, [A, B1, B2, B3])
end

function AB3_A_system(; χN=12.0, ϕAB3=0.5, αAB3=1.0, αA=0.35)
    polymerAB3 = Component(starAB3(), αAB3, ϕAB3)
    polymerA = Component(homopolymer_chain(; label=:A, segment=KuhnSegment(:A)), αA, 1-ϕAB3)
    return PolymerSystem([polymerAB3, polymerA],
                         Dict(Set([:A, :B])=>χN))
end

function model_mp_0()
    system = AB3_A_system()
    uc = UnitCell(4.0)
    lattice = BravaisLattice(uc)
    w = AuxiliaryField(zeros(64), lattice)
    algo = SD(0.2)
    nccscft = NoncyclicChainSCFT(system, w, 0.02; mde=OSF, updater=algo)
    return PolyorderModel(nccscft)
end

function model_mp_1()
    system1 = AB3_A_system()
    uc1 = UnitCell(1.0)
    lattice1 = BravaisLattice(uc1)
    w1 = AuxiliaryField(zeros(1), lattice1)
    algo = SD(0.2)
    nccscft1 = NoncyclicChainSCFT(system1, w1, 0.1; mde=OSF, updater=algo)
    return PolyorderModel(nccscft1)
end

function model_mp_2()
    system2 = AB3_A_system()
    uc2 = UnitCell(4.0)
    lattice2 = BravaisLattice(uc2)
    w2 = AuxiliaryField(zeros(64), lattice2)
    algo = SD(0.2)
    nccscft2 = NoncyclicChainSCFT(system2, w2, 0.02; mde=OSF, updater=algo)
    return PolyorderModel(nccscft2)
end

@testset "macrophase.jl: MacrophaseModel" begin
    ϕ₀ = 0.4
    mpmodel = MacrophaseModel(model_mp_0(), model_mp_0(); ϕ₁=0.2, ϕ₂=0.6, ϕ₀=ϕ₀)

    ϕ₁, ϕ₂, _ = macrophase(JSOOptimizer(), mpmodel; lb=[0.01, ϕ₀], ub=[ϕ₀, 0.99])
    @test ϕ₁ ≈ 0.09284579130809374 atol=1e-4
    @test ϕ₂ ≈ 0.7238098327774354 atol=1e-4

    ϕ₁, ϕ₂, _ = macrophase(OptimOptimizer(), mpmodel; lb=[0.01, ϕ₀], ub=[ϕ₀, 0.99])
    @test ϕ₁ ≈ 0.09285304060762666 atol=1e-4
    @test ϕ₂ ≈ 0.7238078351694467 atol=1e-4
end

@testset "macrophase.jl: SurrogateMacrophaseModel" begin
    ϕ₀ = 0.4
    ϕs_dis1 = [0.01, 0.02, 0.05, 0.1, 0.15]
    ϕs_dis2 = [0.6, 0.65, 0.7, 0.75, 0.8]
    ϕc = ϕControlParameter(1)
    Fs_dis1 = [free_energy(ϕ, ϕc, DISPhase, model_mp_0())[1] for ϕ in ϕs_dis1]
    μs_dis1 = [free_energy(ϕ, ϕc, DISPhase, model_mp_0())[2][1] for ϕ in ϕs_dis1]
    Fs_dis2 = [free_energy(ϕ, ϕc, DISPhase, model_mp_0())[1] for ϕ in ϕs_dis2]
    μs_dis2 = [free_energy(ϕ, ϕc, DISPhase, model_mp_0())[2][1] for ϕ in ϕs_dis2]
    surrogate1 = CHSplineSurrogate(ϕs_dis1, Fs_dis1, μs_dis1)
    surrogate2 = CHSplineSurrogate(ϕs_dis2, Fs_dis2, μs_dis2)
    smpmodel = SurrogateMacrophaseModel(surrogate1, surrogate2, ϕc; ϕ₁=0.1, ϕ₂=0.7, ϕ₀=ϕ₀)

    ϕ₁, ϕ₂, _ = macrophase(OptimOptimizer(), smpmodel; lb=[0.01, ϕ₀], ub=[ϕ₀, 0.99])
    @test ϕ₁ ≈ 0.09200998000796767 atol=1e-4
    @test ϕ₂ ≈ 0.7238331646714469 atol=1e-4

    ϕ₁, ϕ₂, _ = macrophase(JSOOptimizer(), smpmodel; lb=[0.01, ϕ₀], ub=[ϕ₀, 0.99])
    @test ϕ₁ ≈ 0.09196085624309827 atol=1e-4
    @test ϕ₂ ≈ 0.7238331646714469 atol=1e-4

    ϕ₁, ϕ₂, _ = macrophase(RootsOptimizer(), smpmodel)
    @test ϕ₁ ≈ 0.0919736142806544 atol=1e-4
    @test ϕ₂ ≈ 0.7238353843739808 atol=1e-4
end

@testset "macrophase.jl: OPSOptimizer" begin
    vecϕ01 = [0.01, 0.05, 0.25]  # list of initial guesses for phase 1
    vecϕ02 = [0.6, 0.7, 0.8]  # list of initial guesses for phase 2
    tol_io = 1e-8  # tol for internal optimizer
    tol_oo = 1e-4  # tol for outside optimizer
    # ϕ1 = 0.1  # initial guess of the solution for phase 1
    # ϕ2 = 0.7  # initial guess of the solution for pahse 2
    # ϕ0 = 0.4  # the initial volume fractions of the mixing system, critical for optimizing algorithms, but not used by root finding algorithms.

    m = MacrophaseModel(model_mp_1(), model_mp_2(); cached=true)
    io = RootsOptimizer(; x_abstol=tol_io)
    oo = OPSOptimizer(vecϕ01, vecϕ02; optimizer=io, x_abstol=tol_oo)
    ϕ₁, ϕ₂, _ = macrophase(oo, m)
    @test ϕ₁ ≈ 0.09281738092371884 atol=1e-4
    @test ϕ₂ ≈ 0.7238331646714469 atol=1e-4

    m = MacrophaseModel(model_mp_1(), model_mp_2(); cached=true)
    io = OptimOptimizer(; x_abstol=tol_io)
    oo = OPSOptimizer(vecϕ01, vecϕ02; optimizer=io, x_abstol=tol_oo)
    ϕ₁, ϕ₂, _ = macrophase(oo, m)
    @test ϕ₁ ≈ 0.09281738092371884 atol=1e-4
    @test ϕ₂ ≈ 0.7238331646714469 atol=1e-4

    m = MacrophaseModel(model_mp_1(), model_mp_2(); cached=true)
    io = JSOOptimizer(; x_abstol=tol_io)
    oo = OPSOptimizer(vecϕ01, vecϕ02; optimizer=io, x_abstol=tol_oo)
    ϕ₁, ϕ₂, trace = macrophase(oo, m)
    @test ϕ₁ ≈ 0.09281738092371884 atol=1e-4
    @test ϕ₂ ≈ 0.7238331646714469 atol=1e-4
end

nothing