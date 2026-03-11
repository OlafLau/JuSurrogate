module GroundTruth

using Roots
using ..Thermodynamics

export get_ground_truth

function get_ground_truth(α_val::Float64, χN_val::Float64)
    # 寻找 spinodal 点
    spinodals = find_zeros(ϕ -> μ_prime(ϕ, α_val, χN_val), 1e-4, 1.0 - 1e-4)
    if length(spinodals) < 2
        error("未能找到足够的 spinodal 点！")
    end
    ϕ_sp1, ϕ_sp2 = spinodals[1], spinodals[2]

    μ_max = μ(ϕ_sp1, α_val, χN_val)
    μ_min = μ(ϕ_sp2, α_val, χN_val)

    function exact_delta_Fg(m)
        ϕ1 = find_zero(ϕ -> μ(ϕ, α_val, χN_val) - m, (1e-6, ϕ_sp1), Bisection())
        ϕ2 = find_zero(ϕ -> μ(ϕ, α_val, χN_val) - m, (ϕ_sp2, 1.0 - 1e-6), Bisection())
        Fg1 = F_c(ϕ1, α_val, χN_val) - m * ϕ1
        Fg2 = F_c(ϕ2, α_val, χN_val) - m * ϕ2
        return Fg1 - Fg2
    end

    μ_exact = find_zero(exact_delta_Fg, (μ_min + 1e-5, μ_max - 1e-5), Bisection())
    ϕ_alpha = find_zero(ϕ -> μ(ϕ, α_val, χN_val) - μ_exact, (1e-6, ϕ_sp1), Bisection())
    ϕ_beta = find_zero(ϕ -> μ(ϕ, α_val, χN_val) - μ_exact, (ϕ_sp2, 1.0 - 1e-6), Bisection())
    Fg_exact = F_c(ϕ_alpha, α_val, χN_val) - μ_exact * ϕ_alpha

    return (μ=μ_exact, ϕ_α=ϕ_alpha, ϕ_β=ϕ_beta, Fg=Fg_exact, ϕ_sp1=ϕ_sp1, ϕ_sp2=ϕ_sp2)
end

end
