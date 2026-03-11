@testset "phase.jl" begin
    @test DISPhase.description == description(DISPhase)
    @test LAMPhase.variable_name == as_variable_name(LAMPhase)
    @test HEXPhase.ascii_label == as_ascii_label(HEXPhase)
    @test DISPhase.plot_label == as_plot_label(DISPhase)
end