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
nccscft = NoncyclicChainSCFT(system, w, 0.01; mde=OSF)
model = PolyorderModel(nccscft)

## IMPORTANT: model has been modified after each call of free_energy.
# For free_energy with MacrophaseModel and MicrophaseModel, see test_phasemodel.jl

@testset "model.jl: free_enenrgy for ϕ" begin
    system = Polymer.AB_A_system()  # ϕAB = 0.8

    ϕc = ϕControlParameter(1)
    F_dis, μ_dis = free_energy(0.6, ϕc, DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=20.0, ϕAB=0.6, fA=0.5, α=0.5))[1] atol=1e-8

    @test Polymer.ϕ(:AB, system) == 0.5  # system is not modified.
    @test Polymer.ϕ(:AB, nccscft.system) == 0.6  # nccscft has been modified.
    @test getparam(model, ϕc) == 0.6
end

@testset "model.jl: free_enenrgy for α, f, χN" begin
    αc = αControlParameter(2, αParam)
    F_dis, μ_dis = free_energy(0.4, αc, DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=20.0, ϕAB=0.6, fA=0.5, α=0.4))[1] atol=1e-8

    χNc = χNControlParameter(:A, :B, χNParam)
    F_dis, μ_dis = free_energy(18.0, χNc, DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=18.0, ϕAB=0.6, fA=0.5, α=0.4))[1] atol=1e-8

    fc = fControlParameter(1, 1, fParam)
    F_dis, μ_dis = free_energy(0.2, fc, DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=18.0, ϕAB=0.6, fA=0.2, α=0.4))[1] atol=1e-8
end

@testset "free_energy.jl: grand_free_energy" begin
    ϕc = ϕControlParameter(1)
    F_dis, μ_dis = free_energy(0.7, ϕc, DISPhase, model)
    @test F_dis ≈ DIS.F(AB_A_system(χN=18.0, ϕAB=0.7, fA=0.2, α=0.4))[1] atol=1e-8

    @test grand_free_energy(F_dis, μ_dis, [0.7]) ≈ DIS.Fg(polymer_system(model))
end

nothing