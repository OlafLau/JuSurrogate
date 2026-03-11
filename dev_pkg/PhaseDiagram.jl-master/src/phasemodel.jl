"""
AbstractPhaseModel

Provide a unified API interface for storing the infomation for a macrophase or microphase calculation.

All concrete subtypes of AbstractModel except subtypes of SurrogatePhaseModel should support following behaviors:

* `AbstractControlParameter = control_parameter(::AbstractPhaseModel)`
* `AbstractModel = best_model(::AbstractPhaseModel, phase_index, param, paramType)`
* `cache_model!(::AbstractPhaseModel, phase_index, ::AbstractModel)`
"""
abstract type AbstractPhaseModel end
abstract type SurrogatePhaseModel <: AbstractPhaseModel end

"""
    MacrophaseModel <: AbstractPhaseModel

Holds the inital guess (`ϕ₁, ϕ₂`), initial volume fraction (`ϕ₀`), target phases, and models for a macrophase separation computation.

Note that the order subscripts (1 and 2) are important, which stands for the two kinds of phases.
"""
mutable struct MacrophaseModel{T<:Real, R<:ϕControlParameter, P1<:AbstractPhase, P2<:AbstractPhase, M1<:AbstractModel, M2<:AbstractModel} <: AbstractPhaseModel
    ϕ₁::T
    ϕ₂::T
    ϕ₀::T
    ctrlparam::R
    phase1::P1
    phase2::P2
    model1::M1
    model2::M2
    model1_history::Vector{M1}
    model2_history::Vector{M2}
    cached::Bool
    kwargs
end

"""
The default ctrlparam is only for binary system.
"""
function MacrophaseModel(model1::M1, model2::M2,
                         ctrlparam=ϕControlParameter(1);
                         ϕ₁=0.0, ϕ₂=1.0, ϕ₀=0.5*(ϕ₁+ϕ₂),
                         phase1=DISPhase, phase2=DISPhase,
                         cached=false, cache_init=false,
                         kwargs...) where {M1 <: AbstractModel, M2 <: AbstractModel}
    history1 = cache_init ? [clone(model1)] : M1[]
    history2 = cache_init ? [clone(model2)] : M2[]
    return MacrophaseModel(ϕ₁, ϕ₂, ϕ₀, ctrlparam,
                           phase1, phase2, model1, model2,
                           history1, history2,
                           cached, kwargs)
end

control_parameter(mpmodel::MacrophaseModel) = mpmodel.ctrlparam

"""
    SurrogateMacrophaseModel <: SurrogatePhaseModel <: AbstractPhaseModel

Holds the inital guess (`ϕ₁, ϕ₂`), initial volume fraction (`ϕ₀`), target phases, and surrogate functions for computing F and μ.

Note that the order subscripts (1 and 2) are important, which stands for the two kinds of phases.
"""
struct SurrogateMacrophaseModel{T<:Real, R<:ϕControlParameter, P1<:AbstractPhase, P2<:AbstractPhase, S1<:AbstractSurrogate, S2<:AbstractSurrogate} <: SurrogatePhaseModel
    ϕ₁::T
    ϕ₂::T
    ϕ₀::T
    ctrlparam::R
    phase1::P1
    phase2::P2
    surrogate1::S1
    surrogate2::S2
    kwargs
end

function SurrogateMacrophaseModel(surrogate1::S1, surrogate2::S2, ctrlparam=ϕControlParameter(1); ϕ₁=0.0, ϕ₂=1.0, ϕ₀=0.5*(ϕ₁+ϕ₂), phase1=DISPhase, phase2=DISPhase, kwargs...) where {S1<:AbstractSurrogate, S2<:AbstractSurrogate}
    return SurrogateMacrophaseModel(ϕ₁, ϕ₂, ϕ₀, ctrlparam, phase1, phase2, surrogate1, surrogate2, kwargs)
end

"""
    MicrophaseModel <: AbstractPhaseModel

Holds the inital guess (`x₀`), initial bounds `(a, b)`, target phases, and  models for a macrophase separation computation.

`vartype` indicates the type of the independent variable. The canonical free energy is a function of this variable. In other words, we are interested at the microphase separation controlled by this variable. Typical parameters are `fParam` and `χNParam` for single-component polymer system. And `ϕParam` is often used for multi-component polymer system.

Note that the order subscripts (1 and 2) are important, which stands for the two kinds of phases.
"""
struct MicrophaseModel{T<:Real, R<:AbstractControlParameter, P1<:AbstractPhase, P2<:AbstractPhase, M1<:AbstractModel, M2<:AbstractModel} <: AbstractPhaseModel
    a::T
    b::T
    x₀::T
    ctrlparam::R
    phase1::P1
    phase2::P2
    model1::M1
    model2::M2
    model1_history::Vector{M1}
    model2_history::Vector{M2}
    cached::Bool
    kwargs
end

function MicrophaseModel(model1::M1, model2::M2,
                         ctrlparam::AbstractControlParameter;
                         a=0.0, b=1.0, x₀=0.5*(a+b),
                         phase1=DISPhase, phase2=LAMPhase,
                         cached=false, cache_init=false,
                         kwargs...) where {M1<:AbstractModel, M2<:AbstractModel}
    history1 = cache_init ? [clone(model1)] : M1[]
    history2 = cache_init ? [clone(model2)] : M2[]
    return MicrophaseModel(a, b, x₀, ctrlparam,
                           phase1, phase2, model1, model2,
                           history1, history2,
                           cached, kwargs)
end

control_parameter(mpmodel::MicrophaseModel) = mpmodel.ctrlparam

"""
Find a model in the history whose control parameter is closest to `x`.
"""
function best_model(mpmodel::AbstractPhaseModel, phase_index, x)
    cp = control_parameter(mpmodel)
    model = phase_index == 1 ? mpmodel.model1 : mpmodel.model2
    history = phase_index == 1 ? mpmodel.model1_history : mpmodel.model2_history
    mpmodel.cached || return model
    !isempty(history) || return clone(model)

    dist = Inf
    id_best = 0
    for i in eachindex(history)
        new_dist = abs(getparam(history[i], cp) - x)
        if new_dist < dist
            id_best = i
            dist = new_dist
        end
    end

    return clone(history[id_best])
end

function cache_model!(mpmodel::AbstractPhaseModel, phase_index, model::AbstractModel)
    mpmodel.cached || return mpmodel
    history = phase_index == 1 ? mpmodel.model1_history : mpmodel.model2_history

    # when control parameter exists in the history, do not cache again!
    cp = control_parameter(mpmodel)
    xh = getparam(model, cp)
    for m in history
        (getparam(m, cp) ≈ xh) && return mpmodel
    end

    push!(history, clone(model))
    return mpmodel
end

"""
    SurrogateMicrophaseModel <: SurrogatePhaseModel <: AbstractPhaseModel

Holds target phases, configrations, and surrogate functions for computing F and μ.

Note that the order subscripts (1 and 2) are important.
"""
struct SurrogateMicrophaseModel{R<:AbstractControlParameter, P1<:AbstractPhase, P2<:AbstractPhase, S1<:AbstractSurrogate, S2<:AbstractSurrogate} <: SurrogatePhaseModel
    ctrlparam::R
    phase1::P1
    phase2::P2
    surrogate1::S1
    surrogate2::S2
    kwargs
end

function SurrogateMicrophaseModel(surrogate1::S1, surrogate2::S2, ctrlparam::AbstractControlParameter; phase1=DISPhase, phase2=LAMPhase, kwargs...) where {S1<:AbstractSurrogate, S2<:AbstractSurrogate}
    return SurrogateMicrophaseModel(ctrlparam, phase1, phase2, surrogate1, surrogate2, kwargs)
end

# traits for phase separation types
abstract type PhaseSeparationType end
struct MicrophaseSeparation <: PhaseSeparationType end
struct MacrophaseSeparation <: PhaseSeparationType end

phase_separation_type(::AbstractPhaseModel) = MicrophaseSeparation()
phase_separation_type(::MacrophaseModel) = MacrophaseSeparation()
phase_separation_type(::SurrogateMacrophaseModel) = MacrophaseSeparation()
phase_separation_type(::MicrophaseModel) = MicrophaseSeparation()
phase_separation_type(::SurrogateMicrophaseModel) = MicrophaseSeparation()