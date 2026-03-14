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
using Roots
using Optim

params = (αA=1.0, αB=0.5, χN=5.0)

# 1. 精准解析的 Ground Truth 作为误差参照点
true_pa, true_pb, _, _ = solve_gce_root(
    x -> fh_free_energy(x, params.αA, params.αB, params.χN),
    x -> fh_chemical_potential(x, params.αA, params.αB, params.χN),
    collect(range(1e-4, 1.0 - 1e-4, length=50000));
    method=Roots.Bisection()
)
println("==> Exact Ground Truth Binodals: pa=$true_pa, pb=$true_pb\n")

# 2. 选取骨架初始化点
sp_points = fh_spinodals(params.αA, params.αB, params.χN)
sp1, sp2 = sp_points[1], sp_points[2]
p1 = range(1e-4, sp1, length=4)[1:end-1]
p2 = range(sp2, 1.0 - 1e-4, length=4)[2:end]
smart_initial_ϕs = vcat(collect(p1), collect(p2))

# 3. 待考察的插值方法
interpolators = [
    ("Julia_CHS_Surrogate", julia_chs_surrogate),
    # ("Julia_NatNeigh_Surrogate", julia_natural_neighbor_surrogate),
    ("Scipy_Cubic_Surrogate", scipy_cubic_surrogate),
    ("Scipy_Kriging_Surrogate", scipy_kriging_surrogate),
    ("Scipy_IDW_Surrogate", scipy_idw_surrogate),
    # ("Scipy_NatNeigh_Surrogate", scipy_natural_neighbor_surrogate)
]

# 4. 待考察的求解算法与求解途径
solvers = [
    ("Bisection", :gce, Roots.Bisection()),
    ("root(Order1)", :gce, Roots.Order1()),
    ("A42", :gce, Roots.A42()),
    ("Brent", :gce, Roots.Brent()),
    ("BFGS", :ge, Optim.BFGS()),
    ("LBFGS", :ge, Optim.LBFGS()),
    ("GradientDescent", :ge, Optim.GradientDescent()),
    ("ConjugateGradient", :ge, Optim.ConjugateGradient())
]

csv_file = joinpath(@__DIR__, "solvers_benchmark.csv")
open(csv_file, "w") do io
    write(io, "Interpolator,Solver,PredictorType,Status,BlackboxCalls,Residuals\n")

    for (intp_name, intp_func) in interpolators
        for (solv_name, pred_type, solv_algo) in solvers
            println("Testing: Interpolator = [$intp_name], Solver = [$solv_name], Predictor = [$(pred_type)]")

            try
                ans_pred, _, hist = run_surrogate_loop(
                    params;
                    initial_ϕs=smart_initial_ϕs,
                    predictor=pred_type,
                    sampling_strategy="OPS",
                    tol=1e-8,
                    max_iters=100,
                    ground_true=true,
                    verbose=false,
                    solver_algo=solv_algo,
                    surrogate_func=intp_func
                )

                # 提取历史的真实真理残差序列 (使用 Log10)
                res_array = log10.(abs.(hist["phi_a"] .- true_pa) .+ abs.(hist["phi_b"] .- true_pb) .+ 1e-16)

                # 使用引号转义序列里的逗号，避免破坏 CSV
                calls_str = "\"" * string(hist["blackbox_calls"]) * "\""
                res_str = "\"" * string(res_array) * "\""
                status = "Success"

                write(io, "$intp_name,$solv_name,$pred_type,$status,$calls_str,$res_str\n")
                println("  -- [Success] Final Iter: $(hist["iters"]), BB_Calls: $(hist["blackbox_calls"][end])")
            catch e
                println("  -- [Failed]: ", typeof(e))
                status = "Failed"
                write(io, "$intp_name,$solv_name,$pred_type,$status,\"[]\",\"[]\"\n")
            end
        end
    end
end

println("\n[Complete] All algorithm tests finished. Benchmark saved to: ", csv_file)
