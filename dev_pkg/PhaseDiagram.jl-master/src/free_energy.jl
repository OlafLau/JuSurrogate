"""
    free_energy(x, phase_index, mpmodel::AbstractPhaseModel)

Evaluate the free energy and the chemical potential of `mpmmodel.model1` if `phase_index == 1` otherwise `mpmmodel.model2` at control parameter `mpmmodel.ctrlparam` of value `x`.

The model in the model history list which has the closest `x` value will be used to initialized the evaluation. The optimized model will be cached in the model history.
"""
function free_energy(x, phase_index, mpmodel::AbstractPhaseModel)
    phase = phase_index == 1 ? mpmodel.phase1 : mpmodel.phase2
    # the returned model is a deep copy of the element in the history vector.
    model = best_model(mpmodel, phase_index, x)
    # model has been modified here!
    res = free_energy(x, mpmodel.ctrlparam, phase, model)
    mpmodel.cached && cache_model!(mpmodel, phase_index, model)

    return res  # F, μ̃s
end

"""

For ϕ, α, f, χN, b of AbstractControlParameter type.
"""
function free_energy(x, cp::AbstractControlParameter, p::AbstractPhase, model::AbstractModel)
    @info "======================================================"
    @info ""
    @info "\t $(as_variable_name(cp.param)) = $(round(x, digits=6)) \t phase: $(as_variable_name(p))"
    @info ""
    @info "======================================================"

    # Only update the model when parameters are different.
    if !(x ≈ getparam(model, cp))
        setparam!(model, x, cp)
    end

    # Set the output dir for Polyorder.cell_solve!
    # Example: HEX_phi0.02, LAM_f0.1
    # First run config.io.base_dir is ".", following run config.io.base_dir has been modified to "./HEX_phi0.02". We have to first restore the original base_dir.
    pathc = splitpath(model.config.io.base_dir)
    phase_str = as_ascii_label(p)
    var_str = as_ascii_label(cp.param)
    base_dir = phase_str == split(pathc[end], "_")[1] ? joinpath(pathc[1:end-1]) : model.config.io.base_dir
    dir = phase_str * "_" * var_str * "$(round(x, digits=6))"
    dir = joinpath(base_dir, dir)
    mkpath(dir)
    config = model.config
    config = @set config.io.base_dir = dir
    model.config = config

    return free_energy(p, model)  # F, μ̃s
end

function free_energy(ϕ, s::CHSplineSurrogate)
    return s(ϕ), s(ϕ; grad=true)
end

"""
`SplineSurrogate` has no information about μ. We simply return 0.0 to be compatible with `CHSplineSurrogate` version of `free_energy`.
"""
function free_energy(x, s::SplineSurrogate)
    return s(x), 0.0
end

function gibbs_free_energy(ϕ₁, ϕ₂, ϕ₀, F₁, μ₁, F₂, μ₂)
    d = (ϕ₁ - ϕ₂)
    v₁ = (ϕ₀ - ϕ₂) / d
    v₂ = 1 - v₁

    G = v₁ * F₁ + v₂ * F₂

    dG1 = v₁ * (F₂ - F₁) / d + v₁ * μ₁
    dG2 = v₂ * (F₂ - F₁) / d + v₂ * μ₂

    return G, dG1, dG2
end

"""
    gibbs_free_energy

Compute the total Gibbs free energy of a polymer system.

Only 2-specie 2-component polymer system is supported. The computation is performed based on the results of Canonical ensemble computations.
"""
function gibbs_free_energy(ϕ₁, ϕ₂, mpmodel::MacrophaseModel)
    F₁, μ₁ = free_energy(ϕ₁, 1, mpmodel)
    F₂, μ₂ = free_energy(ϕ₂, 2, mpmodel)
    return gibbs_free_energy(ϕ₁, ϕ₂, mpmodel.ϕ₀, F₁, μ₁[1], F₂, μ₂[1])
end

function gibbs_free_energy(ϕ₁, ϕ₂, m::SurrogateMacrophaseModel)
    F₁, μ₁ = free_energy(ϕ₁, m.surrogate1)
    F₂, μ₂ = free_energy(ϕ₂, m.surrogate2)
    return gibbs_free_energy(ϕ₁, ϕ₂, m.ϕ₀, F₁, μ₁[1], F₂, μ₂[1])
end

"""
    grand_free_energy

Compute the grand canonical free energy of a polymer system from canonical free energy and chemical potentials.

Note that both lengths of μs and ϕs are equal to #components - 1. The last component of the PolymerSystem instance is treated implicitly. See Polyorder.jl/chemical_potentials.jl for more details.
"""
function grand_free_energy(F::Real, μs::AbstractVector, ϕs::AbstractVector)
    return F - sum(ϕs .* μs)
end

function grand_free_energy(ϕ, s::CHSplineSurrogate)
    F, μ = free_energy(ϕ, s)
    return grand_free_energy(F, [μ], [ϕ])
end