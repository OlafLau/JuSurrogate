## Roots.Bisection(): 15 pair of evaluations.
# x = 0.02	 phase: HEX
# x = 0.02	 phase: LAM
# x = 0.3	 phase: HEX
# x = 0.3	 phase: LAM
# x = 0.0775	 phase: HEX
# x = 0.0775	 phase: LAM
# x = 0.1525	 phase: HEX
# x = 0.1525	 phase: LAM
# x = 0.21375	 phase: HEX
# x = 0.21375	 phase: LAM
# x = 0.18312499999999998	 phase: HEX
# x = 0.18312499999999998	 phase: LAM
# x = 0.16781249999999998	 phase: HEX
# x = 0.16781249999999998	 phase: LAM
# x = 0.17546874999999998	 phase: HEX
# x = 0.17546874999999998	 phase: LAM
# x = 0.17929687499999997	 phase: HEX
# x = 0.17929687499999997	 phase: LAM
# x = 0.18121093749999997	 phase: HEX
# x = 0.18121093749999997	 phase: LAM
# x = 0.18216796874999996	 phase: HEX
# x = 0.18216796874999996	 phase: LAM
# x = 0.18264648437499997	 phase: HEX
# x = 0.18264648437499997	 phase: LAM
# x = 0.18288574218749998	 phase: HEX
# x = 0.18288574218749998	 phase: LAM
# x = 0.18288574218749995	 phase: HEX
# x = 0.18288574218749995	 phase: LAM
# x = 0.1828857421875	 phase: HEX
# x = 0.1828857421875	 phase: LAM

## Roots.A42(): 8 pair of evaluations
# x = 0.02	 phase: HEX
# x = 0.02	 phase: LAM
# x = 0.3	 phase: HEX
# x = 0.3	 phase: LAM
# x = 0.0775	 phase: HEX
# x = 0.0775	 phase: LAM
# x = 0.1916974006150434	 phase: HEX
# x = 0.1916974006150434	 phase: LAM
# x = 0.1829067169961487	 phase: HEX
# x = 0.1829067169961487	 phase: LAM
# x = 0.1832345163314208	 phase: HEX
# x = 0.1832345163314208	 phase: LAM
# x = 0.18290671699614866	 phase: HEX
# x = 0.18290671699614866	 phase: LAM
# x = 0.18290671699614872	 phase: HEX
# x = 0.18290671699614872	 phase: LAM

## Roots.Brent(): 9 evaluations
# x = 0.02	 phase: HEX
# x = 0.02	 phase: LAM
# x = 0.3	 phase: HEX
# x = 0.3	 phase: LAM
# x = 0.196137466633745	 phase: HEX
# x = 0.196137466633745	 phase: LAM
# x = 0.18175246396380104	 phase: HEX
# x = 0.18175246396380104	 phase: LAM
# x = 0.1830854252661706	 phase: HEX
# x = 0.1830854252661706	 phase: LAM
# x = 0.18306934322490226	 phase: HEX
# x = 0.18306934322490226	 phase: LAM
# x = 0.18306930664700527	 phase: HEX
# x = 0.18306930664700527	 phase: LAM
# x = 0.18306930664700524	 phase: HEX
# x = 0.18306930664700524	 phase: LAM
# x = 0.1830693066470053	 phase: HEX
# x = 0.1830693066470053	 phase: LAM
@testset "microphase.jl: MicrophaseModel" begin
    model = MicrophaseModel("HEX.yml";
                            configfile2="LAM.yml",
                            a=0.02, b=0.3,
                            phase1=HEXPhase, phase2=LAMPhase)

    ϕ = microphase(RootsOptimizer(), model)
    @test ϕ ≈ 0.1828857421875 atol=1e-3
end

@testset "microphase.jl: SurrogateMicrophaseModel" begin
    config_hex = load_config("HEX.yml")
    config_lam = load_config("LAM.yml")
    ϕs_hex, Fs_hex, μs_hex = readfile("free_energy_HEX.txt")
    ϕs_lam, Fs_lam, μs_lam = readfile("free_energy_LAM.txt")
    surrogate_hex = CHSplineSurrogate(ϕs_hex, Fs_hex, μs_hex, config_hex)
    surrogate_lam = CHSplineSurrogate(ϕs_lam, Fs_lam, μs_lam, config_lam)
    model = SurrogateMicrophaseModel(config_hex, surrogate_hex, surrogate_lam;
                                     config2=config_lam,
                                     phase1=HEXPhase, phase2=LAMPhase)

    ϕ = microphase(RootsOptimizer(), model)
    @test ϕ ≈ 0.1831 atol=1e-3
end