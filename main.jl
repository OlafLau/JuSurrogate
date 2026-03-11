include("src/Thermodynamics.jl")
include("src/GroundTruth.jl")
include("src/Surrogates.jl")
include("src/Predictors.jl")
include("src/Sampling.jl")

using .Thermodynamics
using .GroundTruth
using .Surrogates
using .Predictors
using .Sampling

using DataFrames
using Printf

function run_test()
    # 物理参数定义
    const α_val = 0.55
    const χN_val = 5.0

    println("==================================================")
    println("正在计算精确的基准解析点...")
    exact_res = get_ground_truth(α_val, χN_val)

    println("          Ground Truth         ")
    @printf("μ*     = %.8f\n", exact_res.μ)
    @printf("ϕ_α*   = %.8f\n", exact_res.ϕ_α)
    @printf("ϕ_β*   = %.8f\n", exact_res.ϕ_β)
    @printf("Fg*    = %.8f\n", exact_res.Fg)
    println("==================================================\n")

    setup_surrogates()

    models = ["CHS", "CubicSpline", "QuadraticSpline", "Linear"]
    max_iter = 30
    init_ϕ_data = collect(range(0.01, 0.99, length=6))

    df_results = DataFrame()

    for model_name in models
        println(">>> 开始测试代理模型: $model_name")

        ϕ_data = copy(init_ϕ_data)
        F_data = [F_c(p, α_val, χN_val) for p in ϕ_data]
        μ_data = [μ(p, α_val, χN_val) for p in ϕ_data]

        Nf = length(ϕ_data)

        for iter in 1:max_iter
            perm = sortperm(ϕ_data)
            x = ϕ_data[perm]
            y = F_data[perm]
            dydx = μ_data[perm]

            surr_F, surr_mu, surr_Fg = build_surrogate(model_name, x, y, dydx)

            x_range = range(minimum(x), maximum(x), length=1000)
            pred_res = predict_phase_equilibrium(surr_F, surr_mu, x_range)

            error_Fg = isnan(pred_res.Fg) ? NaN : abs(pred_res.Fg - exact_res.Fg)

            push!(df_results, (
                    Model=model_name,
                    Iter=iter,
                    Nf=Nf,
                    Error_Fg=error_Fg,
                    phi_alpha_pred=pred_res.ϕ_α,
                    phi_beta_pred=pred_res.ϕ_β
                ), cols=:union)

            # 收敛标准
            if !isnan(error_Fg) && error_Fg < 1e-8
                println("模型 $model_name 在第 $iter 次迭代时达到收敛精度。")
                break
            end

            new_points = sample_next_points(pred_res.ϕ_α, pred_res.ϕ_β, x, strategy="OPS")

            # 由于可能加入不止一个点（比如OPS会加入两相的值），这里循环添加到数据集中
            for np in new_points
                push!(ϕ_data, np)
                push!(F_data, F_c(np, α_val, χN_val))
                push!(μ_data, μ(np, α_val, χN_val))
                Nf += 1
            end

        end
    end

    # 打印每个模型的结果总结
    for m in models
        println("\n=======================================================")
        println("          模型: $m (OPS策略)       ")
        println("=======================================================")
        sub_df = df_results[df_results.Model.==m, :]
        show(stdout, MIME("text/plain"), sub_df)
        println("\n")
    end
end

run_test()
