function _fg!(FF, GG, ϕ₁, ϕ₂, ϕ₀, χN, f, α, C)
    fout, df1, df2 = G(ϕ₁, ϕ₂, ϕ₀, χN, f, α, C)
    if GG !== nothing
        GG[1], GG[2] = df1, df2
    end
    if FF !== nothing
        return fout
    end
end

"""
Constraint is necessary for most of alogrithms.
SAMIN performs worse than other gradient based algorithms.
Since we only have two variables, BFGS has same performance as L-BFGS.
BFGS is the best algorithm for computing coexistence line of DIS-DIS.
We have to test scipy.optimize.TNC.
There is no TNC implementation in Julia.
"""
function binodal(χN, f, α, C; algo=nothing, opt=nothing, constraint=true)
    ϕc, χNc = critical_point(f, α)
    ϕs1, ϕs2 = spinodal(χNParameter(), χN, f, α)
    x0 = [ϕs1, ϕs2]
    if algo === nothing
        algo = BFGS(; linesearch=Optim.LineSearches.BackTracking())
    end
    if opt === nothing
        opt = Optim.Options(x_tol=1e-4,
                            f_tol=1e-6,
                            f_abstol=1e-6,
                            g_tol=1e-6,
                            g_reltol=1e-6,
                            outer_x_tol=1e-4,
                            outer_f_tol=1e-6,
                            outer_f_abstol=1e-6,
                            outer_f_reltol=1e-6,
                            outer_g_tol=1e-6,
                            outer_g_reltol=1e-6)
    end
    func = Optim.only_fg!((FF, GG, x)->_fg!(FF, GG, x[1], x[2], ϕc, χN, f, α, C))
    if constraint
        lb = [0.000001, ϕc]
        ub = [ϕc, 0.999999]
        if algo isa SAMIN
            result = optimize(func, lb, ub, x0, algo, opt)
        else
            result = optimize(func, lb, ub, x0, Fminbox(algo), opt)
        end
    else
        result = optimize(func, x0, algo, opt)
    end
    ϕb1, ϕb2 = Optim.minimizer(result)
    return ϕb1, ϕb2
end

function _f(x, ϕc, χN, f, α, C)
    out, d1, d2 = G(x[1], x[2], ϕc, χN, f, α, C)
    return out
end

"""
BayesianOptimization takes too many function evaluations to achieve tol <= 1e-4 for x.
We should not use this type of optimization algorithm for computing coexistence point.
"""
function binodal_bopt(χN, f, α, C)
    ϕc, χNc = critical_point(f, α)
    ϕs1, ϕs2 = spinodal(χNParameter(), χN, f, α)
    x0 = [ϕs1, ϕs2]
    lb = [0.000001, ϕc]
    ub = [ϕc, 0.999999]

    func = x -> _f(x, ϕc, χN, f, α, C)

    # opt = BOpt(x -> _f(x, ϕc, χN, f, α, C),
    #         ElasticGPE(2, mean = MeanConst(-10.), kernel = SEArd([0., 0.], 5.),
    #                   logNoise = -2., capacity = 3000),
    #         ExpectedImprovement(),
    #         MAPGPOptimizer(every = 50, noisebounds = [-4, 3],
    #                      kernbounds = [[-1, -1, 0], [4, 4, 10]],
    #                      maxeval = 40),
    #         lb, ub, maxiterations = 1000,
    #        sense = Min, verbosity = Progress)

    # Choose as a model an elastic GP with input dimensions 2.
    # The GP is called elastic, because data can be appended efficiently.
    model = ElasticGPE(2,                            # 2 input dimensions
                    mean = MeanConst(0.),
                    kernel = SEArd([0., 0.], 5.),
                    logNoise = 0.,
                    capacity = 3000)              # the initial capacity of the GP is 3000 samples.
    set_priors!(model.mean, [Normal(1, 2)])

    # Optimize the hyperparameters of the GP using maximum a posteriori (MAP) estimates every 50 steps
    modeloptimizer = MAPGPOptimizer(every = 50, noisebounds = [-4, 3],       # bounds of the logNoise
                                    kernbounds = [[-1, -1, 0], [4, 4, 10]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
                                    maxeval = 40)

    opt = BOpt(func,
            model,
            UpperConfidenceBound(),                   # type of acquisition
            modeloptimizer,
            lb, ub,                     # lowerbounds, upperbounds
            repetitions = 5,                          # evaluate the function for each input 5 times
            maxiterations = 100,                      # evaluate at 100 input positions
            sense = Min,                              # minimize the function
            acquisitionoptions = (method = :LD_LBFGS, # run optimization of acquisition function with NLopts :LD_LBFGS method
                                    restarts = 5,       # run the NLopt method from 5 random initial conditions each time.
                                    maxtime = 0.1,      # run the NLopt method for at most 0.1 second each time
                                    maxeval = 1000),    # run the NLopt methods for at most 1000 iterations (for other options see https://github.com/JuliaOpt/NLopt.jl)
            verbosity = Progress)

    result = boptimize!(opt)
    return result
end

"""
IntervalOptimisation is not working for tol <= 1e-3.
We should not use this type of optimization algorithm for computing coexistence point.
"""
function binodal_interval(χN, f, α, C)
    ϕc, χNc = critical_point(f, α)
    ϕs1, ϕs2 = spinodal(χNParameter(), χN, f, α)
    result = minimise(x -> _f(x, ϕc, χN, f, α, C), (1e-4..ϕc) × (ϕc..1-1e-4), tol=1e-3)
    return result
end