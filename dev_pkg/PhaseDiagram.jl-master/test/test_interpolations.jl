### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 26453814-f4f8-11ea-3c37-055415e2d0f4
using Interpolations

# ╔═╡ 8d895780-f4ff-11ea-3dd1-5d29dbf55853
using DataInterpolations

# ╔═╡ c9157aa2-f4ff-11ea-28b7-19762d90646d
using Plots

# ╔═╡ c35f02dc-f503-11ea-22b8-f5fab649667e
using Polymer

# ╔═╡ c85d29ee-f503-11ea-20d4-e3d2cc988295
αParam.value_type(1)

# ╔═╡ 4a03f290-f55c-11ea-231f-5d1c063111fe
f(x) = exp(x)

# ╔═╡ 9a0c0744-f4ff-11ea-22ef-a726b539816e
x = [0.0, 0.4, 1.0, 1.2, 2.4, 3, π]

# ╔═╡ bf266d12-f4ff-11ea-3933-1b19e02f19b4
y = f.(x)

# ╔═╡ c6fad5dc-f4ff-11ea-20c9-273e2955d2cd
spl1 = DataInterpolations.LinearInterpolation(y, x)

# ╔═╡ e5a128b8-f55b-11ea-0476-a1be962d3a26
spl2 = DataInterpolations.QuadraticInterpolation(y, x)

# ╔═╡ f0bc2d56-f55b-11ea-3ef6-eda274a67e2c
spl22 = DataInterpolations.QuadraticSpline(y, x)

# ╔═╡ 0314c7e0-f55c-11ea-114f-c705c58c9e0d
spl3 = DataInterpolations.CubicSpline(y, x)

# ╔═╡ 1d27381a-f500-11ea-03db-436322ff6e5b
begin
	x_interp = range(0.0, π, step=0.002)
	Plots.scatter(x, y)
	plot!(x_interp, f.(x_interp), label="0")
	plot!(x_interp, spl1.(x_interp), label="1")
	plot!(x_interp, spl2.(x_interp), label="2")
	plot!(x_interp, spl22.(x_interp), label="2-2")
	plot!(x_interp, spl3.(x_interp), label="3")
end

# ╔═╡ af5ac22e-f500-11ea-35af-3f266ee02a38
spl.([0.3, 0.8])

# ╔═╡ eea4ac28-f501-11ea-16e6-f7685bca7ca0
spl.u

# ╔═╡ f2175450-f501-11ea-1377-3f448f0439c8
spl.t

# ╔═╡ Cell order:
# ╠═26453814-f4f8-11ea-3c37-055415e2d0f4
# ╠═8d895780-f4ff-11ea-3dd1-5d29dbf55853
# ╠═c9157aa2-f4ff-11ea-28b7-19762d90646d
# ╠═c35f02dc-f503-11ea-22b8-f5fab649667e
# ╠═c85d29ee-f503-11ea-20d4-e3d2cc988295
# ╠═4a03f290-f55c-11ea-231f-5d1c063111fe
# ╠═9a0c0744-f4ff-11ea-22ef-a726b539816e
# ╠═bf266d12-f4ff-11ea-3933-1b19e02f19b4
# ╠═c6fad5dc-f4ff-11ea-20c9-273e2955d2cd
# ╠═e5a128b8-f55b-11ea-0476-a1be962d3a26
# ╠═f0bc2d56-f55b-11ea-3ef6-eda274a67e2c
# ╠═0314c7e0-f55c-11ea-114f-c705c58c9e0d
# ╠═1d27381a-f500-11ea-03db-436322ff6e5b
# ╠═af5ac22e-f500-11ea-35af-3f266ee02a38
# ╠═eea4ac28-f501-11ea-16e6-f7685bca7ca0
# ╠═f2175450-f501-11ea-1377-3f448f0439c8
