"""
    AbstractModel

Provide a unified API interface for computing the free energy and chemical potentials of a polymer system.

All concrete subtypes of AbstractModel should support following behaviors:

* `PolymerSystem = polymer_system(::AbstractModel)`
* `(T, Vector{T}) = free_energy(::AbstractPhase, ::AbstractModel)`
* `(T, Vector{T}) = free_energy(::Phase{:DIS}, ::AbstractModel)`
"""
abstract type AbstractModel end

function polymer_system end
function free_energy end

"""
    PolyorderModel{S<:Polyorder.AbstractSCFT} <: AbstractModel

Implement a concrete subtype of AbstractModel for Polyorder SCFT computation engine.
"""
mutable struct PolyorderModel{S<:Polyorder.AbstractSCFT, O<:Polyorder.AbstractCellOptAlgorithm} <: AbstractModel
    scft::S
    opt::O
    config::Polyorder.Config
    optimized::Bool
end

function PolyorderModel(scft::T; opt=Polyorder.default(VariableCellOpt), config=Polyorder.Config(), optimized=false) where {T <: Polyorder.AbstractSCFT}
    return PolyorderModel(scft, opt, config, optimized)
end

function clone(model::PolyorderModel)
    newscft = Polyorder.clone(model.scft)
    opt = Polyorder.clone(model.opt)
    return PolyorderModel(newscft, opt, model.config, model.optimized)
end

polymer_system(model::PolyorderModel) = model.scft.system

function free_energy(::AbstractPhase, model::PolyorderModel)
    if !model.optimized
        Polyorder.cell_solve!(model.opt, model.scft, model.config)
        model.optimized = true
    else
        Polyorder.solve!(model.scft, model.config)
    end

    return Polyorder.F(model.scft), Polyorder.μ̃s(model.scft)
end

free_energy(::DISType, model::PolyorderModel) = DIS.F(model.scft.system)

function update!(model::PolyorderModel, system::PolymerSystem)
    newscft = Polyorder.reset(model.scft, system)
    model.scft = newscft
    model.optimized = false
    return model
end

function setparam!(model::PolyorderModel, x::Real, cp::AbstractControlParameter)
    system = Polymer.update!(polymer_system(model), x, cp)
    return update!(model, system)
end

function getparam(model::PolyorderModel, cp::AbstractControlParameter)
    return Polymer.getparam(polymer_system(model), cp)
end

###########################
# set and get ϕParam
# Constraints: all ϕ in a polymer system should sum to 1.
# So we can only update a full list of ϕ.
# Exception: only for binary system a single value can be provided.
###########################

"""
If ϕ is a real number, binary system is assumed, the volume fractions of two components are ϕ and 1-ϕ, respectively.
ϕ can be also a vector. In this case, the multicomponent system with length(ϕ) components is supported.
"""
function setparam!(model::PolyorderModel, ϕ, ::ϕType)
    system = Polymer.update!(polymer_system(model), ϕ, ϕParam)
    return update!(model, system)
end

function getparam(model::PolyorderModel, mol, ::ϕType)
    return Polymer.ϕ(mol, polymer_system(model))
end

function getparam_list(model::PolyorderModel, ::ϕType)
    return Polymer.ϕs(polymer_system(model))
end

######################################################################
# To avoid loading C++ Polyorder, below codes are commented out.

# """
#     PolyorderModel{S<:Polyorder.AbstractSCFT} <: AbstractModel

# Implement a concrete subtype of AbstractModel for C++ Polyorder SCFT computation engine through a wrapper Polyorder.jl.
# """
# mutable struct PolyorderModel{T, S} <: AbstractModel
#     scft::T
#     config::C
#     system::PolymerSystem
# end

# function PolyorderModel(configfile::String, system::PolymerSystem)
#     scft, config = Polyorder.init_model(configfile)
#     return PolyorderModel(scft, config, system)
# end

# polymer_system(model::PolyorderModel) = model.system

# free_energy(::AbstractPhase, model::PolyorderModel) = Polyorder.free_energy!(model.scft, model.config)

# free_energy(::Phase{:DIS}, model::PolyorderModel) = DIS.F(model.system)