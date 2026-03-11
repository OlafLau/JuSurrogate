"""
Compute the free energy of a homogeneous phase (DISPhase). We utilize Polyorder package to compute the free energy and chemical potential of a polymer system. Thanks to Polyorder, now noncyclic chain polymers with arbitrary architecture and multicomponent system are supported.
"""
module DIS

using Polymer
using Polyorder

# Note that, the second returned value is a vector even there is only one element in the vector.
F(ps::PolymerSystem) = (Polyorder.F_DIS(ps), Polyorder.μ̃s_DIS(ps))
Fg(ps::PolymerSystem) = Polyorder.Fg_DIS(ps)

function F(ϕ, χN, f, α, C=1.0)
    ϕA = f * (1 - ϕ) + ϕ
    ϕB = (1 - ϕ) * (1 - f)
    F = χN * ϕA * ϕB
    F += (1 - ϕ) * (log(C * (1 - ϕ)) - 1)
    F += ϕ / α * (log(C * ϕ / α) - 1)
    μ = χN * (1 - f) * (1 - 2ϕA) - log(C*(1 - ϕ))
    μ += log(C * ϕ / α) / α
    μ *= α
    return F, μ
end

function dμ(ϕ, χN, f, α)
    dμ = -2 * χN * (1 - f)^2
    dμ += 1 / (1 - ϕ)
    dμ += 1 / (α * ϕ)
    return dμ
end

function Fg(ϕ, χN, f, α, C=1.0)
    Fc, μ = F(ϕ, χN, f, α, C)
    return Fc - μ * ϕ / α
end

function G(ϕ₁, ϕ₂, ϕ₀, χN, f, α, C=1.0)
    v₁ = (ϕ₀ - ϕ₂) / (ϕ₁ - ϕ₂)
    v₂ = 1 - v₁
    F₁, μ₁ = F(ϕ₁, χN, f, α, C)
    F₂, μ₂ = F(ϕ₂, χN, f, α, C)
    Gout = v₁ * F₁ + v₂ * F₂
    dG1 = v₁ * (F₂ - F₁) / (ϕ₁ - ϕ₂) + v₁ * μ₁ / α
    dG2 = v₂ * (F₂ - F₁) / (ϕ₁ - ϕ₂) + v₂ * μ₂ / α
    return Gout, dG1, dG2
end

function critical_point(f, α)
    ϕc = 1 / ( 1 + sqrt(α))
    χNc = 0.5 * ((1 + sqrt(1/α)) / (1 - f))^2
    return ϕc, χNc
end

function spinodal(::PolymerParameter{:χN}, χN, f, α)
    A = 2χN * α * (1 - f)^2
    δ = (A + 1 - α)^2 - 4A
    if δ < 0
        # This indicates that the input χN is below the critical point χNc.
        # Thus no macrophase separation occurs.
        return nothing
    end
    ϕ1 = 0.5 + (1 - α - sqrt(δ)) / (2A)
    ϕ2 = 0.5 + (1 - α + sqrt(δ)) / (2A)
    return ϕ1, ϕ2
end

function spinodal(::PolymerParameter{:ϕ}, ϕ, f, α)
    χN = 0.5 / (1 - f)^2
    χN *= 1 / (1 - ϕ) + 1 / (α * ϕ)
    return χN
end

# include("optimize_dis.jl")

end # module