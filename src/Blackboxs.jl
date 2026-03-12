module Blackboxs

export evaluate_fh_blackbox, fh_free_energy, fh_chemical_potential

"""
    fh_free_energy(ϕA, αA, αB, χN)

计算二组分(1D) Flory-Huggins 模型的单位体积正则混合自由能 Fc。
为了防止数值溢出，对边界浓度做极小量截断保护。
"""
function fh_free_energy(ϕA::Real, αA::Real, αB::Real, χN::Real)
    # 限制 ϕA 的范围，避免 log(0) 带来的 Inf/NaN
    pA = clamp(ϕA, 1e-15, 1.0 - 1e-15)
    pB = 1.0 - pA
    return (pA * log(pA)) / αA + (pB * log(pB)) / αB + χN * pA * pB
end

"""
    fh_chemical_potential(ϕA, αA, αB, χN)

计算二组分(1D) Flory-Huggins 模型的有效交换化学势 μ (即 ∂Fc/∂ϕA)。
"""
function fh_chemical_potential(ϕA::Real, αA::Real, αB::Real, χN::Real)
    pA = clamp(ϕA, 1e-15, 1.0 - 1e-15)
    pB = 1.0 - pA
    return (1.0 + log(pA)) / αA - (1.0 + log(pB)) / αB + χN * (1.0 - 2.0 * pA)
end

"""
    evaluate_fh_blackbox(ϕ_points, params)

替代真实 SCFT 等昂贵计算的 Black-box 评估模块（以解析的二组分 F-H 理论作为仿真黑盒）。

输入:
- `ϕ_points`: 经过 Sampling 策略后给出的全新的体积配比坐标组（1D Array）。
- `params`: 物理参数，形如具有 αA, αB, χN 等属性的结构体或命名的 Tuple。

输出:
- `Fc_vals`: 当前离散散点处计算得到的正则自由能 F_c 数组。
- `mu_vals`: 当前区域的化学势偏导 μ 数组。
- `Fg_vals`: 根据 F_g = F_c - μ * ϕ 计算出的巨正则自由能 F_g 数组。
"""
function evaluate_fh_blackbox(ϕ_points::AbstractVector{<:Real}, params)
    N = length(ϕ_points)
    Fc_vals = zeros(N)
    mu_vals = zeros(N)
    Fg_vals = zeros(N)
    
    αA = params.αA
    αB = params.αB
    χN = params.χN
    
    for i in 1:N
        ϕ = ϕ_points[i]
        
        # 1. 计算出该点对应的自由能密度的数值 F_c
        fc = fh_free_energy(ϕ, αA, αB, χN)
        
        # 2. 计算出该切线的化学势偏导数值 mu
        mu = fh_chemical_potential(ϕ, αA, αB, χN)
        
        # 3. 后续算出巨正则系综自由能 F_g (在统计力学系综映射下)
        fg = fc - mu * ϕ
        
        Fc_vals[i] = fc
        mu_vals[i] = mu
        Fg_vals[i] = fg
    end
    
    return Fc_vals, mu_vals, Fg_vals
end

end # module
