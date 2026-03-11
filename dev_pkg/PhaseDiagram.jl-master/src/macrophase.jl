function _fix_lb!(lb; a=0.0, b=1.0, ϵ=1e-6)
    lb[1] = lb[1] > a ? lb[1] : a + ϵ
    lb[2] = lb[2] < b ? lb[2] : b + ϵ
    return lb
end

function _fix_ub!(ub; a=0.0, b=1.0, ϵ=1e-6)
    ub[1] = ub[1] < a ? ub[1] : a - ϵ
    ub[2] = ub[2] < b ? ub[2] : b - ϵ
    return ub
end

"""
    macrophase(oo::OptimOptimizer, model::AbstractPhaseModel;
                lb=[1e-6, model.ϕ₀+1e-6], ub=[model.ϕ₀-1e-6, 1.0-1e-6])
    macrophase(jo::JSOOptimizer, model::AbstractPhaseModel;
                lb=[1e-6, model.ϕ₀+1e-6], ub=[model.ϕ₀-1e-6, 1.0-1e-6])
"""
function macrophase(oo::OptimOptimizer, model::AbstractPhaseModel;
                    lb=[1e-6, model.ϕ₀+1e-6], ub=[model.ϕ₀-1e-6, 1.0-1e-6])
    _fix_lb!(lb; b=model.ϕ₀)
    _fix_ub!(ub; a=model.ϕ₀)
    objfun = (F, G, x) -> _fg!(F, G, x, model)
    x0 = [model.ϕ₁, model.ϕ₂]
    result = Optim.optimize(Optim.only_fg!(objfun), lb, ub, x0,
                            Optim.Fminbox(oo.algo), oo.options)
    ϕ1, ϕ2 = Optim.minimizer(result)
    tr = Optim.x_trace(result)
    return ϕ1, ϕ2, (2*length(tr), first.(tr), last.(tr))
end

"""
Due to the limitation of JSOSolver.jl, we cannot obtain the trace of x values for the optimization process. Only printted log is available which may be read directly on the screen.
"""
function macrophase(jo::JSOOptimizer, model::AbstractPhaseModel;
                    lb=[1e-6, model.ϕ₀+1e-6], ub=[model.ϕ₀-1e-6, 1.0-1e-6])
    _fix_lb!(lb; b=model.ϕ₀)
    _fix_ub!(ub; a=model.ϕ₀)
    x0 = [model.ϕ₁, model.ϕ₂]
    nlp = NLPModelsModifiers.LBFGSModel(JSOProblemNoHess(x0, lb, ub, model))
    result = tron(nlp; jo.kwargs...)
    ϕ1, ϕ2 = result.solution
    return ϕ1, ϕ2, (2*result.iter, [], [])
end

"""
    macrophase(::RootsOptimizer, model::SurrogateMacrophaseModel)

Perform a macrophase separation calculation using a surrogate model.

Note that direct usage of a MacrophaseModel is impossible because we lack a method to convert μ to corresponding ϕ without a grand canonical calculation. We assumed all SCFT calculations are done in canonical ensemble.
"""
function macrophase(ro::RootsOptimizer, model::SurrogateMacrophaseModel; kwargs...)
    s1 = model.surrogate1
    s2 = model.surrogate2
    μ2ϕ1 = μ -> μ2ϕ(μ, s1)
    μ2ϕ2 = μ -> μ2ϕ(μ, s2)
    Fg1 = μ -> grand_free_energy(μ2ϕ1(μ), s1)
    Fg2 = μ -> grand_free_energy(μ2ϕ2(μ), s2)
    # Find common range for both phases which make sure the surrogates work.
    μa1 = minimum(s1.spl.dydx)
    μa2 = minimum(s2.spl.dydx)
    μa = max(μa1, μa2)
    μb1 = maximum(s1.spl.dydx)
    μb2 = maximum(s2.spl.dydx)
    μb = min(μb1, μb2)
    μ = Roots.find_zero(μ -> Fg1(μ) - Fg2(μ), (μa, μb), ro.algo; ro.kwargs...)
    # We are not interested in the trace for pure surrogate solver.
    # Therefore, nothing is return for trace.
    return μ2ϕ1(μ), μ2ϕ2(μ), nothing
end

"""
    macrophase(os::AbstractSurrogateSolver, model::MacrophaseModel)

Perform a macrophase separation calculation using a surrogate-based solving algorithm.
"""
function macrophase(os::AbstractSurrogateSolver, mpmodel::MacrophaseModel)
    f1 = x -> (ans = free_energy(x, 1, mpmodel); (ans[1], ans[2][1]))
    f2 = x -> (ans = free_energy(x, 2, mpmodel); (ans[1], ans[2][1]))
    ϕ1, ϕ2, trace = macrophase(os, f1, f2, mpmodel.ϕ₁, mpmodel.ϕ₂, mpmodel.ϕ₀)
    return ϕ1, ϕ2, trace
end

"""
    macrophase(os::OPSOptimizer, model::SurrogateMacrophaseModel)

Perform a macrophase separation calculation using a surrogate-based solving algorithm.
"""
function macrophase(os::OPSOptimizer, mpmodel::MacrophaseModel)
    f1 = x -> (ans = free_energy(x, 1, mpmodel); (ans[1], ans[2][1]))
    f2 = x -> (ans = free_energy(x, 2, mpmodel); (ans[1], ans[2][1]))
    x1 = (first(os.x01) + last(os.x01)) / 2
    x2 = (first(os.x02) + last(os.x02)) / 2
    x0 = (maximum(os.x01) + minimum(os.x02)) / 2
    ϕ1, ϕ2, trace = macrophase(os, f1, f2, x1, x2, x0)
    return ϕ1, ϕ2, trace
end

function macrophase(os::OPSOptimizer, f1, f2, x1, x2, x0)
    vecx1 = deepcopy(os.x01)
    vecx2 = deepcopy(os.x02)
    if isnothing(os.f01)
        vecfg1 = f1.(vecx1)
        vecfg2 = f2.(vecx2)
        vecf1, vecg1 = first.(vecfg1), last.(vecfg1)
        vecf2, vecg2 = first.(vecfg2), last.(vecfg2)
    else
        vecf1, vecg1 = deepcopy(os.f01), deepcopy(os.g01)
        vecf2, vecg2 = deepcopy(os.f02), deepcopy(os.g02)
    end
    vecxans1, vecxans2 = copy(vecx1), copy(vecx2)
    vecfans1, vecfans2 = copy(vecf1), copy(vecf2)

    for i in 1:os.max_itr
        surrogate1 = CHSplineSurrogate(vecx1, vecf1, vecg1)
        surrogate2 = CHSplineSurrogate(vecx2, vecf2, vecg2)

        smpmodel = SurrogateMacrophaseModel(surrogate1, surrogate2;
                                            ϕ₁=x1, ϕ₂=x2, ϕ₀=x0)

        x1, x2, _ = macrophase(os.optimizer, smpmodel)

        push!(vecxans1, x1)
        push!(vecxans2, x2)

        fn1, gn1 = f1(x1)
        push!(vecfans1, fn1)
        push!(vecx1, x1)
        push!(vecf1, fn1)
        push!(vecg1, gn1)
        p1 = sortperm(vecx1)
        vecx1, vecf1, vecg1 = vecx1[p1], vecf1[p1], vecg1[p1]

        fn2, gn2 = f2(x2)
        push!(vecfans2, fn2)
        push!(vecx2, x2)
        push!(vecf2, fn2)
        push!(vecg2, gn2)
        p2 = sortperm(vecx2)
        vecx2, vecf2, vecg2 = vecx2[p2], vecf2[p2], vecg2[p2]

        if (i > 1) && (abs(x1 - vecxans1[end-1]) < os.xatol) && (abs(x2 - vecxans2[end-1]) < os.xatol)
            break
        end
    end

    niter = length(vecxans1) - length(os.x01)
    return vecxans1[end], vecxans2[end], (2*niter, vecxans1, vecxans2, vecfans1, vecfans2)
end