# No Pkg lines

include("src/Predictors.jl")
using .Predictors

params = (αA = 1.0, αB = 0.5, χN = 5.0)

# CE free energy (2 component)
Fc(ϕA, ϕB, αA, αB, χN) = ϕA*log(ϕA)/αA + ϕB*log(ϕB)/αB + χN*ϕA*ϕB
Fc_1d(ϕA) = Fc(ϕA, 1-ϕA, params.αA, params.αB, params.χN)

# First derivative (Chemical potential difference)
ϕ_to_μ(ϕA, ϕB, αA, αB, χN) = (1 + log(ϕA))/αA - (1 + log(1-ϕA))/αB + χN*(1-2ϕA)
mu_c_1d(ϕA) = ϕ_to_μ(ϕA, 1-ϕA, params.αA, params.αB, params.χN)

search_range = collect(range(1e-4, 1.0 - 1e-4, length=500))

println("="^50)
println("Testing the GCE (Root-Finding) Predictor")
println("="^50)
pa_gce, pb_gce, mu_gce, Fg_gce = solve_gce_root(Fc_1d, mu_c_1d, search_range)
println("ϕ_solved:  [$pa_gce, $pb_gce]")
println("μ_solved: $mu_gce")
if !isnan(pa_gce)
    Fg_a = Fc_1d(pa_gce) - mu_gce * pa_gce
    Fg_b = Fc_1d(pb_gce) - mu_gce * pb_gce
    println("Fg in α disorder phase:  $Fg_a")
    println("Fg in β disorder phase:  $Fg_b")
    println("difference between 2 Fg:  $(abs(Fg_a - Fg_b))")
end

println("\n"*"="^50)
println("Testing the GE (Optim) Predictor")
println("="^50)
pa_ge, pb_ge, mu_ge, Fg_ge = solve_ge_optim(Fc_1d, mu_c_1d, search_range)
println("ϕ_solved:  [$pa_ge, $pb_ge]")
println("μ_solved: $mu_ge")
if !isnan(pa_ge)
    Fg_a_ge = Fc_1d(pa_ge) - mu_ge * pa_ge
    Fg_b_ge = Fc_1d(pb_ge) - mu_ge * pb_ge
    println("Fg in α disorder phase:  $Fg_a_ge")
    println("Fg in β disorder phase:  $Fg_b_ge")
    println("difference between 2 Fg:  $(abs(Fg_a_ge - Fg_b_ge))")
end
