# tests for PolyorderModel
using Test
using Polymer
using Scattering
using Polyorder
using PhaseDiagram
import PhaseDiagram: DIS

system = Polymer.AB_A_system(χN=20.0, ϕAB=0.8, fA=0.5, α=0.5)
uc = UnitCell(4.0)
lattice = BravaisLattice(uc)
w = AuxiliaryField(zeros(64), lattice)
nccscft = NoncyclicChainSCFT(system, w, 0.01; λs=[0.1, 0.1, 0.5], mde=OSF)
model = PolyorderModel(nccscft)

@testset "model.jl: free_enenrgy" begin
    F_dis, μ_dis = free_energy(DISPhase, model)
    @test F_dis ≈ DIS.F(system)[1] atol=1e-8
end

@testset "model.jl: getparam, setparam! for ϕ" begin
    @test getparam(model, :AB, ϕParam) == 0.8
    setparam!(model, 0.6, ϕParam)
    @test getparam(model, :AB, ϕParam) == 0.6

    F_dis, μ_dis = free_energy(DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=20.0, ϕAB=0.6, fA=0.5, α=0.5))[1] atol=1e-8
end

@testset "model.jl: getparam, setparam! for α, f, χN, b" begin
    αc = αControlParameter(2, αParam)
    @test getparam(model, αc) == 0.5
    setparam!(model, 0.4, αc)
    @test getparam(model, αc) == 0.4

    χNc = χNControlParameter(:A, :B, χNParam)
    @test getparam(model, χNc) == 20.0
    setparam!(model, 18.0, χNc)
    @test getparam(model, χNc) == 18.0

    fc = fControlParameter(1, 1, fParam)
    @test getparam(model, fc) == 0.5
    setparam!(model, 0.2, fc)
    @test getparam(model, fc) == 0.2

    bc = bControlParameter(:A, bParam)
    @test getparam(model, bc) == 1.0
    setparam!(model, 1.2, bc)
    @test getparam(model, bc) == 1.2

    F_dis, μ_dis = free_energy(DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=18.0, ϕAB=0.6, fA=0.2, α=0.4))[1] atol=1e-8
end

nothing