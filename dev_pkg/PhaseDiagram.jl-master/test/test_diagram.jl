@testset "diagram.jl" begin
    bpd1 = BinaryPhaseDiagram1D(ϕParam, ϕParam.value_type[0, 1], [LAMPhase, DISPhase])
    @test bpd1.xaxis == ϕParam

    bpd2 = BinaryPhaseDiagram2D(ϕParam, χNParam, ϕParam.value_type[0.5, 0.5], χNParam.value_type[10, 30], [LAMPhase, DISPhase])
    @test bpd2.yaxis == χNParam

    bpd3 = BinaryPhaseDiagram3D(ϕParam, χNParam, αParam, ϕParam.value_type[0.5, 0.5], χNParam.value_type[10, 30], αParam.value_type[1,1], [LAMPhase, DISPhase])
    @test bpd3.zaxis == αParam

    tpd = TernaryPhaseDiagram(ϕParam, ϕParam, ϕParam, ϕParam.value_type[0.5, 0.5], ϕParam.value_type[0.1, 0.1], ϕParam.value_type[0.1,0.2], [LAMPhase, DISPhase])
    @test tpd.leftaxis == ϕParam
    @test tpd.rightaxis == ϕParam
    @test tpd.bottomaxis == ϕParam
end