### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ d57ca556-d4a8-11ea-05bb-777af4c887a4
using Optim

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

with minimum funciton value

$$G_{min} = 5.190054565153751$$
"

# ╔═╡ 6bf136d0-d4b2-11ea-1678-d919f18cfd32
algo = Optim.BFGS(; linesearch=Optim.LineSearches.BackTracking())

# ╔═╡ 5c762a12-d4b2-11ea-1bbc-e7941721365f
# options = Optim.Options(x_tol=1e-3, outer_x_tol=1e-3, store_trace=true, show_trace=true);
options = Optim.Options(g_tol=1e-4, outer_g_tol=1e-4, store_trace=true, show_trace=true);

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
	Gmin = 5.190054565153751
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

# ╔═╡ e933c612-d4a8-11ea-10fe-314e26136199
function _fg!(FF, GG, ϕ₁, ϕ₂, ϕ₀, χN, f, α, C)
    fout, df1, df2 = G(ϕ₁, ϕ₂, ϕ₀, χN, f, α, C)
    if GG !== nothing
        GG[1], GG[2] = df1, df2
    end
    if FF !== nothing
        return fout
    end
end

# ╔═╡ 866fe2da-d4b8-11ea-2db3-e7d2a3624c63
func = Optim.only_fg!((FF, GG, x)->_fg!(FF, GG, x[1], x[2], ϕ0, χN, f, α, C))

# ╔═╡ ef2ce48e-d4c8-11ea-2e16-c9cccb189356
res = Optim.optimize(func, lb, ub, x0, Optim.Fminbox(algo), options)

# ╔═╡ 0c403ce6-f2a2-11ea-383c-6db83879f067
Optim.minimum(res)

# ╔═╡ 9e4772ac-f29e-11ea-04de-777e643e20d0
Optim.minimizer(res)

# ╔═╡ fc60d49c-f2a2-11ea-3f95-1182fcfe7a98
Optim.minimum(res) - Gmin

# ╔═╡ d4c97700-e1e4-11ea-3b35-dfe08603f3e1
Optim.minimizer(res) - [ϕ1_solution, ϕ2_solution]

# ╔═╡ feb6c2d0-e1e3-11ea-3ada-b5412aa6a1fb
Optim.trace(res)

# ╔═╡ Cell order:
# ╟─03831256-d4b4-11ea-2127-0705363570eb
# ╠═0c403ce6-f2a2-11ea-383c-6db83879f067
# ╠═9e4772ac-f29e-11ea-04de-777e643e20d0
# ╠═fc60d49c-f2a2-11ea-3f95-1182fcfe7a98
# ╠═d4c97700-e1e4-11ea-3b35-dfe08603f3e1
# ╠═feb6c2d0-e1e3-11ea-3ada-b5412aa6a1fb
# ╠═ef2ce48e-d4c8-11ea-2e16-c9cccb189356
# ╠═6bf136d0-d4b2-11ea-1678-d919f18cfd32
# ╠═5c762a12-d4b2-11ea-1bbc-e7941721365f
# ╠═866fe2da-d4b8-11ea-2db3-e7d2a3624c63
# ╠═ab3759ee-d4aa-11ea-1aa8-83f62bae1bf7
# ╠═38e93bb8-d4a9-11ea-13f4-cff7953948e4
# ╠═e933c612-d4a8-11ea-10fe-314e26136199
# ╠═c5084cc2-d4a8-11ea-0bf4-5981f3b82263
# ╠═ad21736e-d4a7-11ea-00d9-016218eff1ee
# ╠═d57ca556-d4a8-11ea-05bb-777af4c887a4
