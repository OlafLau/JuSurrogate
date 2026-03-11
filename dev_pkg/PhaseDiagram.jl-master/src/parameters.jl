"""
    PolyorderConfig <: AbstractConfig
    PolymerConfig <: AbstractConfig

Traits for referring to the types of configuration file.
* `CppPolyorderConfig`: C++ Polyorder configuration file (v11.0 and later).
* `PolymerConfig`: Polymer.jl configuration file.

## Important Notes
* Only 2-specie 2-component polymer system is supported.
* `setparam!` for is not implemented. It involves complicated interactions with the chain architecture.
"""
abstract type AbstractConfig end
struct CppPolyorderConfig <: AbstractConfig end
struct PolymerConfig <: AbstractConfig end
struct PolyorderConfig <: AbstractConfig end

getparam(config, ::PolymerParameter{:ϕ}, ::CppPolyorderConfig) = config["Model"]["phi"][2]

function setparam!(config, ϕ, ::PolymerParameter{:ϕ}, ::PolyorderConfig)
    config["Model"]["phi"][1] = one(ϕ) - ϕ
    config["Model"]["phi"][2] = ϕ
    return config
end

getparam(config, ::PolymerParameter{:χN}, ::PolyorderConfig) = config["Model"]["chiN"][1]

function setparam!(config, χN, ::PolymerParameter{:χN}, ::PolyorderConfig)
    config["Model"]["chiN"][1] = χN
    return config
end

getparam(config, ::PolymerParameter{:f}, ::PolyorderConfig) = config["Model"]["f"][1]

# function setparam!(config, f, ::PolymerParameter{:f}, ::PolyorderConfig)
#     config["Model"]["f"][1] = f
# end

function getparam(config, ::PolymerParameter{:α}, ::PolyorderConfig)
    fs = config["Model"]["f"]
    return fs[1] * fs[end]
end

function setparam!(config, α, ::PolymerParameter{:α}, ::PolyorderConfig)
    f = getparam(config, fParam)
    config["Model"]["f"][end] = α / f
    return config
end

"""
    getparam(config, p::AbstractParameter, configtype=PolyorderConfig())

Get the interested parameter value from a Dict represention of a configuration file whose type is default to `PolyorderConfig`.

Note that, 2-specie 2-chain system is currently implictly assumed. This shall be extended to allow other polymer systems in the future.
"""
getparam(config, p::AbstractParameter, configtype=PolyorderConfig()) = getparam(config, p, configtype)
setparam!(config, val, p::AbstractParameter, configtype=PolyorderConfig()) = setparam!(config, val, p, configtype)