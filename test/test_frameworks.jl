# 设置环境变量或从宿主获取
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/Surrogates.jl")
include("../src/Predictors.jl")
include("../src/Sampling.jl")
include("../src/Blackboxs.jl")
include("../src/Frameworks.jl")

using .Surrogates
using .Predictors
using .Sampling
using .Blackboxs
using .Frameworks
using Plots
using Roots

params = (αA = 1.0, αB = 0.5, χN = 5.0)

# 1. 寻找精准的 Ground Truth 作为误差参照点（高精度的目标真实答案）
true_pa, true_pb, _, _ = solve_gce_root(
    x -> fh_free_energy(x, params.αA, params.αB, params.χN),
    x -> fh_chemical_potential(x, params.αA, params.αB, params.χN),
    collect(range(1e-4, 1.0 - 1e-4, length=50000));
    method=Roots.Bisection()
)
println("==> Exact Ground Truth Binodals: pa=$true_pa, pb=$true_pb")

# 2. 算出行列式的界限 Spinodals 用于打点初始化
sp_points = fh_spinodals(params.αA, params.αB, params.χN)
sp1, sp2 = sp_points[1], sp_points[2]
p1 = range(1e-4, sp1, length=4)[1:end-1]  
p2 = range(sp2, 1.0 - 1e-4, length=4)[2:end] 
smart_initial_ϕs = vcat(collect(p1), collect(p2))

# 3. 分别测试三种抽样策略在相同求根器下的代价表现
println("="^50)
println("Running: GCE + OPS (On-Position Sampling)")
ans_ops, _, hist_ops = run_surrogate_loop(params; initial_ϕs = smart_initial_ϕs, predictor=:gce, sampling_strategy="OPS", verbose=false)

println("Running: GCE + BS (Bisection Sampling)")
ans_bs, _, hist_bs = run_surrogate_loop(params; initial_ϕs = smart_initial_ϕs, predictor=:gce, sampling_strategy="BS", verbose=false)

println("Running: GCE + QS (Quad-section Sampling)")
ans_qs, _, hist_qs = run_surrogate_loop(params; initial_ϕs = smart_initial_ϕs, predictor=:gce, sampling_strategy="QS", verbose=false)

println("Running: GCE + RS (Random Sampling)")
ans_rs, _, hist_rs = run_surrogate_loop(params; initial_ϕs = smart_initial_ϕs, predictor=:gce, sampling_strategy="RS", verbose=false)

# 4. 计算误差残差 Residual
function calc_residual(hist, true_a, true_b)
    res = abs.(hist["phi_a"] .- true_a) .+ abs.(hist["phi_b"] .- true_b)
    # 取 log 之前加上一个机器精度级极小值防止 log(0)
    return log10.(res .+ 1e-16)
end

res_ops = calc_residual(hist_ops, true_pa, true_pb)
res_bs  = calc_residual(hist_bs, true_pa, true_pb)
res_qs  = calc_residual(hist_qs, true_pa, true_pb)
res_rs  = calc_residual(hist_rs, true_pa, true_pb)

# 5. 绘制并保存收敛图
p = plot(title="Convergence Trace Comparison", 
         xlabel="Total Black Box Function Calls", 
         ylabel="Log10(Residual Error)",
         grid=true, legend=:topright,
         dpi=300)

plot!(p, hist_ops["blackbox_calls"], res_ops, label="OPS (On-Position)", marker=:circle, markersize=5, lw=2)
plot!(p, hist_bs["blackbox_calls"], res_bs, label="BS (Bisection)", marker=:square, markersize=5, lw=2)
plot!(p, hist_qs["blackbox_calls"], res_qs, label="QS (Quad-section)", marker=:utriangle, markersize=5, lw=2)
plot!(p, hist_rs["blackbox_calls"], res_rs, label="RS (Random)", marker=:star, markersize=5, lw=2)

save_path = joinpath(@__DIR__, "convergence_trace.png")
savefig(p, save_path)
println("\n[Success] Convergence plot generated and saved to: $save_path")
