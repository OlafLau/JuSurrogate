# tests for PolyorderModel
using Test
using CubicHermiteSpline
using Optim
using Roots
import PhaseDiagram: DIS
using DelimitedFiles
using Random: Xoshiro
using Polyorder
using Polyorder: Polyorder, Anderson, SD, BB, VariableCell

function model_lam()
    ab = AB_system(; χN=15.0, fA=0.5)
    lat = BravaisLattice(UnitCell(3.6))
    # prepare the config file.
    io = IOConfig(base_dir="microphase")
    scftconfig = SCFTConfig(; max_iter=1000)
    config = Polyorder.Config(; io, scft=scftconfig)
    # prepare a well converged SCFT instance via a stable SD method.
    scft0 = NoncyclicChainSCFT(ab, lat, 0.01; init=:randn, rng=Xoshiro(3083))
    Polyorder.solve!(scft0, config)
    # prepare the SCFT instance with an acceleration method.
    updater = Anderson(20; warmup=10)
    scft = NoncyclicChainSCFT(ab, lat, 0.01; updater, init=:randn, rng=Xoshiro(3083))
    Polyorder.initialize!(scft, scft0.wfields)

    opt = Polyorder.default(VariableCellOpt)
    return PolyorderModel(scft, opt, config, false)
end

function model_hex()
    ab = AB_system(; χN=15.0, fA=0.4)
    uc = UnitCell(Hexagonal2D(), 4.0)
    lat = BravaisLattice(uc, 17)
    # prepare the config file.
    io = IOConfig(base_dir="microphase")
    scftconfig = SCFTConfig(; max_iter=500)
    config = Polyorder.Config(; io, scft=scftconfig)
    # prepare a well converged SCFT instance.
    scft0 = NoncyclicChainSCFT(ab, lat, 0.01; init=:randn, rng=Xoshiro(3083))
    Polyorder.solve!(scft0, config)
    # prepare the SCFT instance with acceleration method.
    updater = Anderson(20; warmup=10)
    scft = NoncyclicChainSCFT(ab, lat, 0.01; updater, init=:randn, rng=Xoshiro(3083))
    Polyorder.initialize!(scft, scft0.wfields)

    opt = Polyorder.default(VariableCellOpt)
    return PolyorderModel(scft, opt, config, false)
end

@testset "microphase.jl: OPSOptimizer" begin
    fc = fControlParameter(1, 1, fParam)
    mipm = MicrophaseModel(model_lam(), model_hex(), fc, a=0.36, b=0.48, phase1=LAMPhase, phase2=HEXPhase, cached=true)
    ops = OPSOptimizer([0.36, 0.42, 0.48], [0.36, 0.42, 0.48])
    f, trace = microphase(ops, mipm)
    @test f ≈ 0.39172697337490636 atol=1e-4
end

@testset "microphase.jl: RootsOptimizer" begin
    fc = fControlParameter(1, 1, fParam)
    mipm = MicrophaseModel(model_lam(), model_hex(), fc, a=0.36, b=0.48, phase1=LAMPhase, phase2=HEXPhase, cached=true)
    f, tracker = microphase(RootsOptimizer(), mipm)
    @test f ≈ 0.3917319614414286 atol=1e-4
end

@testset "microphase.jl: SurrogateMicrophaseModel" begin
    data = open("free_energy_microphase.csv", "r") do io
        readdlm(io)
    end
    fs, Fs_lam, Fs_hex = data[:,1], data[:,2], data[:,3]
    s1 = SplineSurrogate(fs, Fs_lam, k=2)
    s2 = SplineSurrogate(fs, Fs_hex, k=2)
    fc = fControlParameter(1, 1, fParam)
    smipm = SurrogateMicrophaseModel(s1, s2, fc; phase1=LAMPhase, phase2=HEXPhase)
    f, _ = microphase(RootsOptimizer(), smipm)
    @test f ≈ 0.3922044925634358 atol=1e-4
    # free_energy_microphase.csv is produced by the following code
    # function model_lam()
    #     system2 = AB_system(; χN=15.0, fA=0.5)
    #     uc2 = UnitCell(3.6)
    #     lattice2 = BravaisLattice(uc2)
    #     w2 = AuxiliaryField(zeros(32), lattice2)
    #     nccscft2 = NoncyclicChainSCFT(system2, w2, 0.01; λs=[0.1, 0.1, 0.5])
    #     config = Polyorder.Config(; cellopt=CellOptConfig(; interval=1.0))
    #     PolyorderModel(nccscft2, config, false)
    # end
    # function model_hex()
    #     system3 = AB_system(; χN=15.0, fA=0.5)
    #     uc3 = UnitCell(HexRect(), 4.2)  # for ϕhA=0.05, a = 4.23
    #     lattice3 = BravaisLattice(uc3)
    #     w3 = AuxiliaryField(zeros(28, 49), lattice3)
    #     nccscft3 = NoncyclicChainSCFT(system3, w3, 0.01; λs=[0.1, 0.1, 0.5])
    #     nccscft3.wfields[1] .= wA
    #     nccscft3.wfields[2] .= wB
    #     nccscft3.wfields[3] .= η
    #     config = Polyorder.Config(; cellopt=CellOptConfig(; interval=1.0))
    #     PolyorderModel(nccscft3, config, false)
    # end
    # fc = fControlParameter(1, 1, fParam)
    # mipm = MicrophaseModel(model_lam(), model_hex(), fc; phase1=LAMPhase, phase2=HEXPhase, cached=true)
    # fs = [0.36, 0.39, 0.42, 0.45, 0.48]
    # fg_lam = [free_energy(f, 1, mipm) for f in fs]
    # Fs_lam = first.(fg_lam)
    # fg_hex = [free_energy(f, 2, mipm) for f in fs]
    # Fs_hex = first.(fg_hex)
end

nothing