@testset "surrogate.jl" begin
    ϕs_hex, Fs_hex, μs_hex = collect(0.1:0.1:1), rand(10), rand(10)
    chs_hex = CHSplineSurrogate(ϕs_hex, Fs_hex, μs_hex)
    @test chs_hex(0.1) ≈ chs_hex.y[1]
    @test chs_hex(0.2) ≈ chs_hex.y[2]
    @test chs_hex.x[1] == 0.1
    @test chs_hex.x[end] == 1
    @test chs_hex.spl isa CubicHermiteSplineInterpolation

    ss_hex = SplineSurrogate(ϕs_hex, Fs_hex)
    @test ss_hex(0.1) ≈ ss_hex.y[1]
    @test ss_hex(0.2) ≈ ss_hex.y[2]
    @test ss_hex.x[1] == 0.1
    @test ss_hex.x[end] == 1
end

nothing
