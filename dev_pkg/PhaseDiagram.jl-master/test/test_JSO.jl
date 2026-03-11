### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ d57ca556-d4a8-11ea-05bb-777af4c887a4
using NLPModels, JSOSolvers

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

$$\chi N=19$$

$$f=0.4$$

$$\alpha=0.14$$

$$C=1.0$$

A suitable initial condition is

$$\phi_0 = 0.7$$

And the minimizer of the function $$G$$ is

$$\phi_1 = 0.631055$$

$$\phi_2 = 0.812722$$

with minimum function value

$$G_{min} = 5.190054565153751$$
"

# ╔═╡ ab3759ee-d4aa-11ea-1aa8-83f62bae1bf7
begin
	ϕ0 = 0.7
  	x0 = [0.5, 0.9]
	lb = [0.000001, ϕ0]
  	ub = [ϕ0, 0.999999]
end

# ╔═╡ 38e93bb8-d4a9-11ea-13f4-cff7953948e4
begin
	χN = 19
	f = 0.4
	α = 0.14
	C = 1.0
	# the accuracy of the following solution is about 1e-4
	ϕ1_solution = 0.631055
	ϕ2_solution = 0.812722
end

# ╔═╡ ad21736e-d4a7-11ea-00d9-016218eff1ee
function F(ϕ, χN, f, α, C)
	println("------------ F IS EVALUATED ONCE ------------")
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

# ╔═╡ c21acaee-e1e6-11ea-1bbe-bf2da5d7543f
begin
	
struct MyProblemNoHess <: AbstractNLPModel
  meta::NLPModelMeta
  counters::Counters
end

function MyProblemNoHess()
  meta = NLPModelMeta(
    2, # nvar
    x0=x0,
    lvar=lb,
    uvar=ub
  )
  MyProblemNoHess(meta, Counters())
end

function NLPModels.obj(nlp::MyProblemNoHess, x::AbstractVector)
  v, _, _ = G(x[1], x[2], ϕ0, χN, f, α, C)
  return v
end

function NLPModels.grad!(nlp::MyProblemNoHess, x::AbstractVector, g::AbstractVector)
  _, g[1], g[2] = G(x[1], x[2], ϕ0, χN, f, α, C)
  return g
end
	
end

# ╔═╡ 866fe2da-d4b8-11ea-2db3-e7d2a3624c63
nlp = LBFGSModel(MyProblemNoHess())

# ╔═╡ 6bf136d0-d4b2-11ea-1678-d919f18cfd32
res = tron(nlp; atol=1e-4)

# ╔═╡ fa18d91a-f2a1-11ea-1b58-25f5dfb3d2b9
res.objective

# ╔═╡ 3e8d9aca-f29f-11ea-01ec-19edeb019835
res.solution

# ╔═╡ 493fc4ae-f2a2-11ea-19cb-f586bfe278d9
res.solution - [ϕ1_solution, ϕ2_solution]

# ╔═╡ Cell order:
# ╟─03831256-d4b4-11ea-2127-0705363570eb
# ╠═fa18d91a-f2a1-11ea-1b58-25f5dfb3d2b9
# ╠═3e8d9aca-f29f-11ea-01ec-19edeb019835
# ╠═493fc4ae-f2a2-11ea-19cb-f586bfe278d9
# ╠═6bf136d0-d4b2-11ea-1678-d919f18cfd32
# ╠═866fe2da-d4b8-11ea-2db3-e7d2a3624c63
# ╠═c21acaee-e1e6-11ea-1bbe-bf2da5d7543f
# ╠═ab3759ee-d4aa-11ea-1aa8-83f62bae1bf7
# ╠═38e93bb8-d4a9-11ea-13f4-cff7953948e4
# ╠═c5084cc2-d4a8-11ea-0bf4-5981f3b82263
# ╠═ad21736e-d4a7-11ea-00d9-016218eff1ee
# ╠═d57ca556-d4a8-11ea-05bb-777af4c887a4
