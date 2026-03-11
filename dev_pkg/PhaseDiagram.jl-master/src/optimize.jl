abstract type AbstractOptimizer end

struct OptimOptimizer <: AbstractOptimizer
    algo::Optim.AbstractOptimizer
    options::Optim.Options
end

function OptimOptimizer(; algo=Optim.BFGS(; linesearch=Optim.LineSearches.BackTracking()), x_abstol=1e-4, max_iter=100, kwargs...)
    # Use callback function to terminate the optimization once x_abstol is arrived.
    # Otherwise, some algorithms in Optim will only use g_tol as convergence test and x_abstol is completely ignored.
    cb = tr -> begin
                (tr[end].iteration == 0) && (return false)
                x_prev = tr[end-1].metadata["x"]
                x = tr[end].metadata["x"]
                converged = (sum(abs.(x .- x_prev) .< x_abstol) == length(x))
               end
    # Enforce store_trace and extended_trace to be true.
    kwargs = (; kwargs..., store_trace=true, extended_trace=true, iterations=max_iter)
    options = Optim.Options(; callback=cb, kwargs...)

    return OptimOptimizer(algo, options)
end

struct JSOOptimizer <: AbstractOptimizer
    kwargs
end

function JSOOptimizer(; x_abstol=1e-4, max_iter=100, kwargs...)
    # Enforce use_only_objgrad to be true.
    kwargs = (; kwargs..., atol=x_abstol, max_eval=max_iter, use_only_objgrad=true)
    return JSOOptimizer(kwargs)
end

struct RootsOptimizer <: AbstractOptimizer
    algo
    kwargs
end

"""
Possible Algorithms:
* ITP: recommended
* A42
* Brent
* Bisection

Possible kwargs for Roots.jl are:
* :xatol, :xabstol
* :xrtol, :xreltol
* :atol,  :abstol
* :rtol,  :reltol
* :maxevals,   :maxsteps
* :strict
"""
function RootsOptimizer(; algo=Roots.ITP(), x_abstol=1e-3, max_iter=20, kwargs...)
    kwargs =(; kwargs..., xatol=x_abstol, maxevals=max_iter)
    return RootsOptimizer(algo, kwargs)
end

## Tools for Optim

function _fg!(F, G, x, model::AbstractPhaseModel)
    f, df1, df2 = gibbs_free_energy(x[1], x[2], model)
    if G !== nothing
        G[1], G[2] = df1, df2
    end
    if F !== nothing
        return f
    end
end

## Tools for JSOSolvers

struct JSOProblemNoHess{T, S} <: AbstractNLPModel{T, S}
    meta::NLPModelMeta{T, S}
    counters::Counters
    model::AbstractPhaseModel
end

function JSOProblemNoHess(x0, lb, ub, model::AbstractPhaseModel)
    meta = NLPModelMeta(2, x0=x0, lvar=lb, uvar=ub)
    return JSOProblemNoHess(meta, Counters(), model)
end

function NLPModels.obj(nlp::JSOProblemNoHess, x::AbstractVector)
    f, _, _ = gibbs_free_energy(x[1], x[2], nlp.model)
    return f
end

function NLPModels.grad!(nlp::JSOProblemNoHess, x::AbstractVector, g::AbstractVector)
    _, g[1], g[2] = gibbs_free_energy(x[1], x[2], nlp.model)
    return g
end

function NLPModels.objgrad!(nlp::JSOProblemNoHess, x::AbstractVector, g::AbstractVector)
    f, g[1], g[2] = gibbs_free_energy(x[1], x[2], nlp.model)
    return f, g
end

## Tools for RootsOptimizer
function ϕ2μ(ϕ, s::CHSplineSurrogate)
    return free_energy(ϕ, s)[2][1]
end

"""
Note that spl.graident(x) is assumed to be a mono-decrease or mono-increase function. Otherwise, there are two or more solutions for this function and may lead to the failing of bracketing the solution by simply assume the brackets of (spl.x[1], spl.x[end]).
"""
function μ2ϕ(μ, s::CHSplineSurrogate)
    if s.spl.dydx[1] == μ
        return s.spl.x[1]
    end
    if s.spl.dydx[end] == μ
        return s.spl.x[end]
    end
    return Roots.find_zero(ϕ -> ϕ2μ(ϕ, s) - μ, (s.spl.x[1], s.spl.x[end]), Roots.A42())
end

"""
AbstractSurrogateSolver optimizers based on the surrogate of the free energy function.
"""
abstract type AbstractSurrogateSolver <: AbstractOptimizer end

"""
IMPORTANT Assumptions:

* the solutions: x1 < x2.
* Both x01 and x02 is sorted and in the ascending order.
* The range of x01 should contain x1.
* The range of x02 should contain x2.
* f01, f02, g01, g02 are optional.
"""
struct OPSOptimizer{O<:AbstractOptimizer} <: AbstractSurrogateSolver
    optimizer::O
    max_itr::Int64
    xatol::Float64
    x01
    x02
    f01
    f02
    g01
    g02
    # make sure that x01, x02 are mandated and its size should be >= 2.
    function OPSOptimizer(x01, x02; optimizer=RootsOptimizer(), max_itr=50, x_abstol=1e-4, f01=nothing, f02=nothing, g01=nothing, g02=nothing)
        length(x01) >= 2 || error("Number of elements in x01 must >= 2!")
        length(x02) >= 2 || error("Number of elements in x02 must >= 2!")

        return new{typeof(optimizer)}(optimizer, max_itr, x_abstol, x01, x02, f01, f02, g01, g02)
    end
end
