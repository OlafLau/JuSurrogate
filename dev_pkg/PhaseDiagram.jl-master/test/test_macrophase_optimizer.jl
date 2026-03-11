### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 373f2f84-f49a-11ea-17b8-1b09c349ddf7
using DelimitedFiles

# ╔═╡ 00681ccc-f49b-11ea-194a-1b11a9f20caf
using Roots

# ╔═╡ 394352a4-f49c-11ea-1555-c7b6c45092c4
using LaTeXStrings

# ╔═╡ 09690460-f49b-11ea-0f67-4d0cd4ad4c79
using Plots

# ╔═╡ 5e45f79c-f49b-11ea-1cf9-839c25eea263
using CubicHermiteSpline

# ╔═╡ ac75ea78-f49b-11ea-3c81-51e0a486772d
using Polymer

# ╔═╡ 82826ad2-f49b-11ea-0b48-f9853d0aefb6
using PhaseDiagram

# ╔═╡ 21b5f688-f49b-11ea-083e-2f6537159475
function readfile(infile)
    data = readdlm(infile, '\t', Float64, '\n'; comments=true)
  # phih  F    mu
    return data[:, 1], data[:, 2], data[:, 3]
end

# ╔═╡ 0502cdf8-f49e-11ea-0d11-7334a38dd1b9
begin
	Fg(F, μ, ϕ, α) = F - μ * ϕ / α
	
	function Fg(ϕ, s::CHSplineSurrogate, config)
		F, μ = free_energy(ϕ, s, config)
		α = getparam(config, αParam)
		return Fg(F, μ, ϕ, α)
	end
end

# ╔═╡ b37cdbdc-f4b3-11ea-09a7-a1d320d99c17
function ϕ2μ(ϕ, s::CHSplineSurrogate, config)
	return free_energy(ϕ, s, config)[2]
end

# ╔═╡ e328857e-f4a5-11ea-3fde-d9c115e0a8fb
function μ2ϕ(μ, s::CHSplineSurrogate, config)
	α = getparam(config, αParam)
	if s.spl.gradient[1] == μ/α
		return s.spl.x[1]
	end
	if s.spl.gradient[end] == μ/α
		return s.spl.x[end]
	end
	return Roots.find_zero(ϕ -> ϕ2μ(ϕ, s, config) - μ,
							(s.spl.x[1], s.spl.x[end]),
							Roots.Bisection())
end

# ╔═╡ 7caa36aa-f4b2-11ea-3de6-5d707accea41
function Fgμ(μ, s::CHSplineSurrogate, config)
	return Fg(μ2ϕ(μ, s, config), s, config)
end

# ╔═╡ efb0be98-f4a6-11ea-1a8d-553be6ee31ed
function test_μ2ϕ(ϕ, surrogate, config)
	μ = ϕ2μ(ϕ, surrogate, config)
	return μ2ϕ(μ, surrogate, config) - ϕ
end

# ╔═╡ 4b776542-f49b-11ea-169b-17a96a8aa5c4
begin
	config_hex = load_config("HEX.yml")
    config_lam = load_config("LAM.yml")
    ϕs_hex, Fs_hex, μs_hex = readfile("free_energy_HEX.txt")
    ϕs_lam, Fs_lam, μs_lam = readfile("free_energy_LAM.txt")
    surrogate_hex = CHSplineSurrogate(ϕs_hex, Fs_hex, μs_hex, config_hex)
    surrogate_lam = CHSplineSurrogate(ϕs_lam, Fs_lam, μs_lam, config_lam)
end

# ╔═╡ 87670440-f4a8-11ea-29a6-0103a6845a22
test_μ2ϕ(0.188, surrogate_hex, config_hex)

# ╔═╡ 107dcfc2-f4a9-11ea-19a1-bf9cb7ad78e2
test_μ2ϕ(0.188, surrogate_lam, config_lam)

# ╔═╡ a0e13abe-f4b2-11ea-0312-ed746a569c4d
begin
	Fgdiff = μ -> Fgμ(μ, surrogate_hex, config_hex) - Fgμ(μ, surrogate_lam, config_lam)
	μa = max(minimum(μs_hex), minimum(μs_lam))
	μb = min(maximum(μs_hex), maximum(μs_lam))
	μ_solution = Roots.find_zero(Fgdiff, (μa, μb), Roots.Bisection())
	ϕ_solution_hex = μ2ϕ(μ_solution, surrogate_hex, config_hex)
	ϕ_solution_lam = μ2ϕ(μ_solution, surrogate_lam, config_lam)
	μ_solution, ϕ_solution_hex, ϕ_solution_lam
end

# ╔═╡ a8b9e892-f49b-11ea-0826-0b87a1fea73d
begin
	ϕ_interp_hex = range(ϕs_hex[1], ϕs_hex[end], step=0.002)
	F_interp_hex = [free_energy(ϕ, surrogate_hex, config_hex)[1] for ϕ in ϕ_interp_hex]
	Plots.scatter(ϕs_hex, Fs_hex;
				xlabel=L"\phi", ylabel=L"F", label="HEX", legend=:topleft)
	plot!(ϕ_interp_hex, F_interp_hex, label="interpolated")
	
	ϕ_interp_lam = range(ϕs_lam[1], ϕs_lam[end], step=0.002)
	F_interp_lam = [free_energy(ϕ, surrogate_lam, config_lam)[1] for ϕ in ϕ_interp_lam]
	Plots.scatter!(ϕs_lam, Fs_lam;
				xlabel=L"\phi", ylabel=L"F", label="LAM", legend=:topleft)
	plot!(ϕ_interp_lam, F_interp_lam, label="interpolated")
	Plots.vline!([ϕ_solution_hex])
	Plots.vline!([ϕ_solution_lam])
end

# ╔═╡ 57e1ad3c-f49c-11ea-2881-e1633b045432
begin
	μ_interp_hex = [free_energy(ϕ, surrogate_hex, config_hex)[2] for ϕ in ϕ_interp_hex]
	Plots.scatter(ϕs_hex, μs_hex;
				xlabel=L"\phi", ylabel=L"\mu", label="HEX", legend=:bottomright)
	plot!(ϕ_interp_hex, μ_interp_hex, label="interpolated")
	
	μ_interp_lam = [free_energy(ϕ, surrogate_lam, config_lam)[2] for ϕ in ϕ_interp_lam]
	Plots.scatter!(ϕs_lam, μs_lam;
				xlabel=L"\phi", ylabel=L"\mu", label="LAM", legend=:bottomright)
	plot!(ϕ_interp_lam, μ_interp_lam, label="interpolated")
	Plots.hline!([μ_solution])
	Plots.vline!([ϕ_solution_hex])
	Plots.vline!([ϕ_solution_lam])
end

# ╔═╡ d892fa3e-f49d-11ea-150f-254220c94ee4
begin
	α = getparam(config_hex, αParam)
	Fgs_hex = Fs_hex .- μs_hex .* ϕs_hex / α
	Fgs_lam = Fs_lam .- μs_lam .* ϕs_lam / α
	Fg_interp_hex = [Fg(ϕ, surrogate_hex, config_hex) for ϕ in ϕ_interp_hex]
	Fg_interp_lam = [Fg(ϕ, surrogate_lam, config_lam) for ϕ in ϕ_interp_lam]
	
	Plots.scatter(μs_hex, Fgs_hex;
				xlabel=L"\mu", ylabel=L"F_g", label="HEX", legend=:bottomleft)
	plot!(μ_interp_hex, Fg_interp_hex, label="interpolated")
	Plots.scatter!(μs_lam, Fgs_lam;
				xlabel=L"\mu", ylabel=L"F_g", label="LAM", legend=:bottomleft)
	plot!(μ_interp_lam, Fg_interp_lam, label="interpolated")
	Plots.vline!([μ_solution])
end

# ╔═╡ e2a1d634-f56b-11ea-2717-c72fb9e42b06


# ╔═╡ Cell order:
# ╠═373f2f84-f49a-11ea-17b8-1b09c349ddf7
# ╠═00681ccc-f49b-11ea-194a-1b11a9f20caf
# ╠═394352a4-f49c-11ea-1555-c7b6c45092c4
# ╠═09690460-f49b-11ea-0f67-4d0cd4ad4c79
# ╠═5e45f79c-f49b-11ea-1cf9-839c25eea263
# ╠═ac75ea78-f49b-11ea-3c81-51e0a486772d
# ╠═82826ad2-f49b-11ea-0b48-f9853d0aefb6
# ╠═21b5f688-f49b-11ea-083e-2f6537159475
# ╠═0502cdf8-f49e-11ea-0d11-7334a38dd1b9
# ╠═7caa36aa-f4b2-11ea-3de6-5d707accea41
# ╠═b37cdbdc-f4b3-11ea-09a7-a1d320d99c17
# ╠═e328857e-f4a5-11ea-3fde-d9c115e0a8fb
# ╠═efb0be98-f4a6-11ea-1a8d-553be6ee31ed
# ╠═87670440-f4a8-11ea-29a6-0103a6845a22
# ╠═107dcfc2-f4a9-11ea-19a1-bf9cb7ad78e2
# ╠═4b776542-f49b-11ea-169b-17a96a8aa5c4
# ╠═a8b9e892-f49b-11ea-0826-0b87a1fea73d
# ╠═57e1ad3c-f49c-11ea-2881-e1633b045432
# ╠═d892fa3e-f49d-11ea-150f-254220c94ee4
# ╠═a0e13abe-f4b2-11ea-0312-ed746a569c4d
# ╠═e2a1d634-f56b-11ea-2717-c72fb9e42b06
