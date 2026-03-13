module Frameworks

export run_surrogate_loop

using ..Surrogates
using ..Predictors
using ..Sampling
using ..Blackboxs

"""
    run_surrogate_loop(params; kwargs...)

完整的代理模型加速自洽场（或黑箱）相图寻根循环框架。

## 参数
- `params`: 黑箱函数（如 F-H 解析模型或真实 SCFT）所需的参数（例如 (αA=1.0, αB=0.5, χN=5.0)）

## 关键字参数 (kwargs)
- `initial_ϕs`: 初始打点位置（默认为 [0.001, 0.999] 或者均匀粗网格）
- `max_iters`: 最大迭代次数（默认 50）
- `tol`: 两次迭代组分间距的收敛判断容差（默认 1e-6）
- `surrogate_func`: 代理模型生成器，接受 (coords, values; grad=true)，默认为 `julia_chs_surrogate`
- `predictor`: 求解器选择，支持 `:gce` (全域求根) 和 `:ge` (全局极小化)，默认为 `:ge`
- `sampling_strategy`: 采样策略，支持 `"OPS"`, `"BS"`, `"TS"`, `"QS"`, `"RS"`，默认为 `"OPS"`
- `search_range`: Predictor 工作的连续搜索域（默认 range(1e-4, 1.0 - 1e-4, length=1000)）
- `verbose`: 是否打印每步迭代信息（默认 true）

## 返回值
- `(ϕ_pred, Fg_pred, history)`
  - `ϕ_pred`: 最终收敛的 (ϕ_α, ϕ_β) 元组
  - `Fg_pred`: 最终环境的巨正则自由能 Fg
  - `history`: 字典，保存了迭代过程的指标 `history["phi_a"]` 等。
"""
function run_surrogate_loop(params;
    initial_ϕs=[0.01, 0.5, 0.99],
    max_iters=50,
    tol=1e-6,
    surrogate_func=julia_chs_surrogate,
    predictor=:ge,
    sampling_strategy="OPS",
    search_range=collect(range(1e-4, 1.0 - 1e-4, length=1000)),
    verbose=true)

    # 历史档案：保存已经被 Blackbox 评估过的真实已知点 (ϕ_i -> Fc, mu)
    ϕ_evaluated = Float64[]
    Fc_evaluated = Float64[]
    mu_evaluated = Float64[]

    # 追踪预测轨迹
    history = Dict{String, Any}(
        "phi_a" => Float64[],
        "phi_b" => Float64[],
        "mu" => Float64[],
        "iters" => 0,
        "blackbox_calls" => Int[]
    )

    # 上一轮的预测点（用于测试收敛状态）
    last_pa, last_pb = NaN, NaN
    pa_pred, pb_pred = NaN, NaN
    mu_pred, Fg_pred = NaN, NaN

    # 第一波：对用户给定的 initial_ϕs 进行真实的黑盒计算初始化
    ϕ_next = copy(initial_ϕs)

    for iter in 1:max_iters
        # 1. 评估 (Black-Box)
        # 将刚刚生成的打点请求发送给真实评估器
        # (过滤掉已经存在于档案的重复打点)
        new_ϕ_to_eval = Float64[]
        for pt in ϕ_next
            if isempty(ϕ_evaluated) || minimum(abs.(ϕ_evaluated .- pt)) > 1e-8
                push!(new_ϕ_to_eval, pt)
            end
        end

        if !isempty(new_ϕ_to_eval)
            Fc_new, mu_new, _ = evaluate_fh_blackbox(new_ϕ_to_eval, params)
            append!(ϕ_evaluated, new_ϕ_to_eval)
            append!(Fc_evaluated, Fc_new)
            append!(mu_evaluated, mu_new)

            # 为了让 1D 代理插值器正常工作，需要对评估域按坐标从小到大重新排序
            perm = sortperm(ϕ_evaluated)
            ϕ_evaluated = ϕ_evaluated[perm]
            Fc_evaluated = Fc_evaluated[perm]
            mu_evaluated = mu_evaluated[perm]
        end

        # 2. 拟合 (Surrogate)
        # 我们把 (Fc, mu) 打包成 tuple 喂给 surrogate_func，这样就能使用梯度的 Hermes 样条等模型
        F_c_surr, mu_c_surr = surrogate_func(ϕ_evaluated, (Fc_evaluated, mu_evaluated); grad=true)

        # 3. 求解 (Prediction)
        if predictor == :gce
            pa_pred, pb_pred, mu_pred, Fg_pred = solve_gce_root(F_c_surr, mu_c_surr, search_range)
        elseif predictor == :ge
            pa_pred, pb_pred, mu_pred, Fg_pred = solve_ge_optim(F_c_surr, mu_c_surr, search_range)
        else
            error("Unknown predictor: $predictor. Use :gce or :ge")
        end

        # 记录踪迹
        push!(history["phi_a"], pa_pred)
        push!(history["phi_b"], pb_pred)
        push!(history["mu"], mu_pred)
        push!(history["blackbox_calls"], length(ϕ_evaluated))
        history["iters"] = iter

        if verbose
            if isnan(pa_pred)
                println("Iter $iter: Predictor failed to find two phases (Single Phase or NaN generated).")
            else
                println("Iter $iter: Predicted ϕ_α=$(round(pa_pred, digits=6)), ϕ_β=$(round(pb_pred, digits=6)), μ=$mu_pred")
            end
        end

        # 4. 判断收敛判定 (Convergence Check)
        if !isnan(pa_pred) && !isnan(pb_pred) && !isnan(last_pa) && !isnan(last_pb)
            if abs(pa_pred - last_pa) < tol && abs(pb_pred - last_pb) < tol
                if verbose
                    println("Converged at iteration $(iter)!")
                end
                break
            end
        end

        last_pa = pa_pred
        last_pb = pb_pred

        # 5. 生成新样本点 (Sampling)
        # 根据预测到的共存点（即使失败了也会触发对应的全局或随机 fallback）生成下一轮请求点
        ϕ_next = sample_next_points([pa_pred, pb_pred], ϕ_evaluated, strategy=sampling_strategy)

        if verbose
            println(" -> Sampling logic asks for Black-box evaluation on: $ϕ_next")
        end
    end

    return (pa_pred, pb_pred), Fg_pred, history
end

end # module
