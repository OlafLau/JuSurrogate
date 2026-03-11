using SnoopPrecompile

@precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    using Logging
    system = AB_A_system(; χN=16.0)
    uc = UnitCell(1.0)
    lattice = BravaisLattice(uc)
    scft = NoncyclicChainSCFT(system, lattice, 0.1)
    fc = fControlParameter(1, 1, fParam)
    ioconfig = IOConfig(; verbosity=-1, save_w=false, save_ϕ=false)
    scftconfig = SCFTConfig(; min_iter=1, max_iter=1, tolmode=:F, tol=1.0)
    cellconfig = CellOptConfig(; max_iter=1)
    config = Polyorder.Config(; io=ioconfig, scft=scftconfig, cellopt=cellconfig)
    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        with_logger(NullLogger()) do
            model = PolyorderModel(scft, config, false)
            mpmodel = MacrophaseModel(model, model; ϕ₁=0.1, ϕ₂=0.7, ϕ₀=0.4)
            ϕ₁, ϕ₂, _ = macrophase(OptimOptimizer(), mpmodel)

            mipm = MicrophaseModel(model, model, fc, a=0.36, b=0.48, phase1=LAMPhase, phase2=HEXPhase, cached=true)
            ro = RootsOptimizer(; x_abstol=0.01, max_iter=1)
            f, _= microphase(ro, mipm)
        end
    end
end