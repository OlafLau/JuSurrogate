### A Pluto.jl notebook ###
# v0.11.2

using Markdown
using InteractiveUtils

# ╔═╡ 03831256-d4b4-11ea-2127-0705363570eb
md"The aim is to optimize the following function

$$G(\phi_1, \phi_2) = v_1(\phi_1,\phi_2) F_1(\phi_1) + v_2(\phi_1,\phi_2) F_2(\phi_2)$$

where

$$v_1(\phi_1, \phi_2) = \frac{\phi_0 - \phi_2}{\phi_1 - \phi_2}$$

$$v_2(\phi_1, \phi_2) = \frac{\phi_1 - \phi_0}{\phi_1 - \phi_2}$$

$$F(\phi) = \chi N\phi_A\phi_B + (1-\phi)\lbrace\ln[C(1-\phi)] - 1\rbrace + \frac{\phi}{\alpha}\left( \ln\frac{C\phi}{\alpha} - 1 \right)$$

$$\phi_A = (1-\phi)f + \phi$$

$$\phi_B = 1 - \phi_A$$

Geometrically, if $$F$$ is a double well function (with two minima and one maximum in between), then the minimum of the function $$G$$ corresponds to the two contact points at the common tangent line of the function $$F$$.

For

$$f=0.0$$

$$\alpha=1.0$$

$$C=1.0$$

A suitable initial condition is

$$\phi_0 = 0.5$$

And the minimizer of the function $$G$$ is

$$\phi_1 = 0.31464218298980107$$

$$\phi_2 = 0.6853578170101975$$
"

# ╔═╡ ef2ce48e-d4c8-11ea-2e16-c9cccb189356


# ╔═╡ 866fe2da-d4b8-11ea-2db3-e7d2a3624c63
md"
tol > 1e-3 gives no useful results while tol=1e-3 can not converge in minutes. Actually I lost my patience to wait for the answer.
"

# ╔═╡ d57ca556-d4a8-11ea-05bb-777af4c887a4
using IntervalArithmetic, IntervalOptimisation

# ╔═╡ e58a749e-d4c8-11ea-13b7-15bee499a024
using Plots

# ╔═╡ 38e93bb8-d4a9-11ea-13f4-cff7953948e4
begin
	χN = 2.1
	f = 0.0
	α = 1.0
	C = 1.0
	ϕ1_solution = 0.31464218298980107
	ϕ2_solution = 0.6853578170101975
end

# ╔═╡ ad21736e-d4a7-11ea-00d9-016218eff1ee
function F(ϕ, χN, f, α, C)
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

# ╔═╡ c5084cc2-d4a8-11ea-0bf4-5981f3b82263
function G(ϕ₁, ϕ₂, ϕ₀, χN, f, α, C)
    v₁ = (ϕ₀ - ϕ₂) / (ϕ₁ - ϕ₂)
    v₂ = 1 - v₁
    F₁, μ₁ = F(ϕ₁, χN, f, α, C)
    F₂, μ₂ = F(ϕ₂, χN, f, α, C)
    Gout = v₁ * F₁ + v₂ * F₂
    dG1 = v₁ * (F₂ - F₁) / (ϕ₁ - ϕ₂) + v₁ * μ₁ / α
    dG2 = v₂ * (F₂ - F₁) / (ϕ₁ - ϕ₂) + v₂ * μ₂ / α
    return Gout, dG1, dG2
end

# ╔═╡ e933c612-d4a8-11ea-10fe-314e26136199
function func(x, ϕc, χN, f, α, C)
    out, d1, d2 = G(x[1], x[2], ϕc, χN, f, α, C)
    return out
end

# ╔═╡ ab3759ee-d4aa-11ea-1aa8-83f62bae1bf7
global_min, minimisers = minimise(x -> func(x, 0.5, χN, f, α, C), (0.1..0.49) × (0.51..0.9), tol=1e-2)

# ╔═╡ 6bf136d0-d4b2-11ea-1678-d919f18cfd32
minimisers

# ╔═╡ 5c762a12-d4b2-11ea-1bbc-e7941721365f
global_min

# ╔═╡ Cell order:
# ╠═03831256-d4b4-11ea-2127-0705363570eb
# ╠═ef2ce48e-d4c8-11ea-2e16-c9cccb189356
# ╠═6bf136d0-d4b2-11ea-1678-d919f18cfd32
# ╠═5c762a12-d4b2-11ea-1bbc-e7941721365f
# ╟─866fe2da-d4b8-11ea-2db3-e7d2a3624c63
# ╠═ab3759ee-d4aa-11ea-1aa8-83f62bae1bf7
# ╠═38e93bb8-d4a9-11ea-13f4-cff7953948e4
# ╠═e933c612-d4a8-11ea-10fe-314e26136199
# ╠═c5084cc2-d4a8-11ea-0bf4-5981f3b82263
# ╠═ad21736e-d4a7-11ea-00d9-016218eff1ee
# ╠═d57ca556-d4a8-11ea-05bb-777af4c887a4
# ╠═e58a749e-d4c8-11ea-13b7-15bee499a024
