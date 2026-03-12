module Predictors

using Roots
using Optim

export solve_gce_root, solve_ge_optim

"""
    solve_gce_root(F_c, mu_c, search_range; method=Bisection())
使用巨正则系综(GCE)条件通过求根解两相平衡。
基于共同有效化学势和相间的巨正则自由能相等 (F_g^α = F_g^β)。
F_c: 代理模型返回的自由能函数
mu_c: 代理模型返回的偏导（化学势）函数
search_range: 搜索范围，比如 range(0.01, 0.99, length=100)。
method: 指定求解根的方法，如 Bisection() (默认), A42(), Brent() 等。
返回: (ϕ_α, ϕ_β, μ, Fg) 或在无解时返回 NaN 元组。
"""
function solve_gce_root(F_c, mu_c, search_range; method=Bisection())
    # 取消依赖于 Spinodal 的特定范围检测，因为对于类似于多个有序抛物线的连续域，它可能存在多个波动
    # 直接在全球域寻找对于任意给定化学势 m 时存在的各个根。
    
    mus = mu_c.(search_range)
    mu_min = minimum(mus)
    mu_max = maximum(mus)
    
    if mu_max <= mu_min || isnan(mu_min)
        return (NaN, NaN, NaN, NaN)
    end
    
    # 我们需要在区间 (mu_min, mu_max) 里找到使得 delta_Fg(m) = 0 的化学势
    m_tests = range(mu_min + 1e-5, mu_max - 1e-5, length=100)

    # 目标函数：让在给定化学势 m 下两端的 巨正则自由能 Fg 相等
    function delta_Fg(m)
        rts = find_zeros(p -> mu_c(p) - m, search_range[1], search_range[end])
        if length(rts) >= 2
            # 第一相和最后一相的最远根
            pa = rts[1]
            pb = rts[end]
            Fg_a = F_c(pa) - m * pa
            Fg_b = F_c(pb) - m * pb
            return Fg_a - Fg_b
        else
            return NaN
        end
    end

    dF_vals = Float64[]
    for mt in m_tests
        try
            push!(dF_vals, delta_Fg(mt))
        catch
            push!(dF_vals, NaN)
        end
    end

    bracket = nothing
    for i in 1:length(dF_vals)-1
        if !isnan(dF_vals[i]) && !isnan(dF_vals[i+1]) && dF_vals[i] * dF_vals[i+1] <= 0
            bracket = (m_tests[i], m_tests[i+1])
            break
        end
    end

    if bracket !== nothing
        try
            mu_pred = find_zero(delta_Fg, bracket, method)
            
            # 使用算出的化学势找回精准体积占比
            rts = find_zeros(p -> mu_c(p) - mu_pred, search_range[1], search_range[end])
            if length(rts) >= 2
                pa_pred = rts[1]
                pb_pred = rts[end]
                Fg_pred = F_c(pa_pred) - mu_pred * pa_pred
                return (pa_pred, pb_pred, mu_pred, Fg_pred)
            end
        catch
        end
    end
    
    return (NaN, NaN, NaN, NaN)
end

"""
    solve_ge_optim(F_c, mu_c, search_range; initial_guess=nothing, method=BFGS())
使用 Gibbs 系综(GE)条件通过 Optim 极小化总自由能(Lever Rule)。
F_c: 代理模型返回的自由能函数
mu_c: 代理模型返回的偏导（化学势）函数
search_range: 提供域边界以设置约束（比如 range(0.01, 0.99, length=100)）
method: 指定内部优化算法，支持 BFGS() (默认), LBFGS(), GradientDescent(), ConjugateGradient() 等。
返回: (ϕ_α, ϕ_β, μ, Fg) 或在无解时返回 NaN 元组。
"""
function solve_ge_optim(F_c, mu_c, search_range; initial_guess=nothing, method=BFGS())
    if initial_guess === nothing
        # 默认初猜
        L = search_range[end] - search_range[1]
        pa0 = search_range[1] + L * 0.1
        pb0 = search_range[end] - L * 0.1
        initial_guess = [pa0, pb0]
    end
    
    # 彻底释放强制绑定于某相空间域（抛弃强制分为左右）的想法
    # 但根据 Gibbs 物理学约束，我们需要有一个混合态支点保证能量最小化求解不变成发散的外推推断
    ϕ_0 = (initial_guess[1] + initial_guess[2]) / 2.0
    
    lower = [search_range[1], search_range[1]]
    upper = [search_range[end], search_range[end]]
    
    # 优化目标：混合系统总体自由能 G
    function G(x)
        pa, pb = x[1], x[2]
        if pa >= ϕ_0 || pb <= ϕ_0 || pa >= pb
            return Inf
        end
        va = (pb - ϕ_0) / (pb - pa)
        vb = 1.0 - va
        return va * F_c(pa) + vb * F_c(pb)
    end
    
    # 使用 mu_c 构建杠杆定律带来的准确目标函数偏导
    function g!(G_grad, x)
        pa, pb = x[1], x[2]
        if pa >= ϕ_0 || pb <= ϕ_0 || pa >= pb
            G_grad[1] = 0.0
            G_grad[2] = 0.0
            return
        end
        va = (pb - ϕ_0) / (pb - pa)
        vb = 1.0 - va
        slope = (F_c(pb) - F_c(pa)) / (pb - pa)
        
        G_grad[1] = va * (mu_c(pa) - slope)
        G_grad[2] = vb * (mu_c(pb) - slope)
    end

    # 保证初猜在安全的盒子约束内，不再强分中位，单纯受限于搜索边界
    initial_guess[1] = clamp(initial_guess[1], lower[1] + 1e-5, ϕ_0 - 1e-5)
    initial_guess[2] = clamp(initial_guess[2], ϕ_0 + 1e-5, upper[2] - 1e-5)

    try
        res = optimize(G, g!, lower, upper, initial_guess, Fminbox(method), Optim.Options(g_abstol=1e-10))
        
        if Optim.converged(res)
            pa_pred, pb_pred = Optim.minimizer(res)
            # 通过连接共切线推算化学势
            mu_pred = (F_c(pb_pred) - F_c(pa_pred)) / (pb_pred - pa_pred)
            Fg_pred = F_c(pa_pred) - mu_pred * pa_pred
            return (pa_pred, pb_pred, mu_pred, Fg_pred)
        end
    catch
    end
    
    return (NaN, NaN, NaN, NaN)
end

end # module
