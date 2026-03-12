using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test
using Optim
using Roots

include(joinpath(@__DIR__, "../src/Predictors.jl"))
using .Predictors

@testset "Predictors Tests" begin
    # 构造一个简单的双波势函数模型 
    # F(x) = (x^2 - 1)^2
    # 极小值在 x = -1 和 x = 1，公共切线在 F = 0 处，化学势 μ = 0
    # 所以相平衡点（Spinodals中间的反曲点）应该是 -1 和 1

    F_c(x) = (x^2 - 1.0)^2
    mu_c(x) = 4.0 * x * (x^2 - 1.0)

    search_range = collect(range(-1.5, 1.5, length=100))

    @testset "GCE Root" begin
        pa, pb, mu, Fg = solve_gce_root(F_c, mu_c, search_range)
        @test isapprox(pa, -1.0, atol=1e-3)
        @test isapprox(pb, 1.0, atol=1e-3)
        @test isapprox(mu, 0.0, atol=1e-3)
        @test isapprox(Fg, 0.0, atol=1e-3)
    end

    @testset "GE Optim" begin
        pa, pb, mu, Fg = solve_ge_optim(F_c, mu_c, search_range)
        @test isapprox(pa, -1.0, atol=1e-3)
        @test isapprox(pb, 1.0, atol=1e-3)
        @test isapprox(mu, 0.0, atol=1e-3)
        @test isapprox(Fg, 0.0, atol=1e-3)
    end

    @testset "No Roots / Single Phase" begin
        # 如果是一个完全单调的向上开口的抛物线如 x^2，就不存在两相共存
        F_c2(x) = x^2
        mu_c2(x) = 2.0x

        pa_r, pb_r, mu_r, Fg_r = solve_gce_root(F_c2, mu_c2, search_range)
        @test isnan(pa_r)

        pa_o, pb_o, mu_o, Fg_o = solve_ge_optim(F_c2, mu_c2, search_range)
        @test isnan(pa_o) || (abs(pa_o - pb_o) < 1e-3) # 优化器可能会陷入同一点，导致无法分离两相
    end
end
