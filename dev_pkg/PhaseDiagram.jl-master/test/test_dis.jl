using PhaseDiagram: DIS

function constants()
    f = 0.4
    α = 0.2
    C = 1.0
    return f, α, C
end

function constants_symmetric_ABblend()
    f = 0
    α = 1
    C = 1
    return f, α, C
end

@testset "dis.jl: free energy" begin
    χN = 20.0
    ϕ = 0.5
    f = 0.5
    α = 0.5
    # AB is treated as independent component
    system = Polymer.AB_A_system(; χN=χN, ϕAB=ϕ, fA=f, α=α)
    F_JP, μ_JP = DIS.F(system)
    # A is treated as independent component
    F, μ = DIS.F(ϕ, χN, f, α)
    @test F_JP == F
    @test μ_JP[1] == -μ / α
end

function binodal_symmetric_ABblend(ϕ)
    return log(ϕ / (1 - ϕ)) / (2ϕ - 1)
end

@testset "dis.jl: critical point" begin
    f, α, C = constants_symmetric_ABblend()
    ϕc, χNc = DIS.critical_point(f, α)
    @test ϕc ≈ 0.5 rtol=1e-12
    @test χNc ≈ 2.0 rtol=1e-12
end

@testset "dis.jl: spinodal" begin
    f, α, C = constants()
    χN = 16
    ϕ1, ϕ2 = DIS.spinodal(χNParam, χN, f, α)
    χN1 = DIS.spinodal(ϕParam, ϕ1, f, α)
    χN2 = DIS.spinodal(ϕParam, ϕ2, f, α)
    @test χN1 ≈ χN rtol=1e-12
    @test χN2 ≈ χN rtol=1e-12
end

# @testset "dis.jl: binodal" begin
#     f, α, C = constants_symmetric_ABblend()
#     χN = 2.1
#     ϕb1, ϕb2 = binodal(χN, f, α, C)
#     @test ϕb1 + ϕb2 ≈ 1.0 rtol=1e-6

#     χN1 = binodal_symmetric_ABblend(ϕb1)
#     χN2 = binodal_symmetric_ABblend(ϕb2)
#     @test χN1 ≈ χN rtol=1e-4
#     @test χN2 ≈ χN rtol=1e-4
# end