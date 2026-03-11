module Thermodynamics

export F_c, μ, μ_prime, F_g

function F_c(ϕ, α, χN)
    ϕ = clamp(ϕ, 1e-12, 1.0 - 1e-12)
    return ϕ * log(ϕ) + (1 - ϕ) * log(1 - ϕ) / α + χN * ϕ * (1 - ϕ)
end

function μ(ϕ, α, χN)
    ϕ = clamp(ϕ, 1e-12, 1.0 - 1e-12)
    return 1 + log(ϕ) - (1 + log(1 - ϕ)) / α + χN * (1 - 2ϕ)
end

function μ_prime(ϕ, α, χN)
    ϕ = clamp(ϕ, 1e-12, 1.0 - 1e-12)
    return 1 / ϕ + 1 / (α * (1 - ϕ)) - 2χN
end

function F_g(ϕ, α, χN)
    return F_c(ϕ, α, χN) - μ(ϕ, α, χN) * ϕ
end

end
