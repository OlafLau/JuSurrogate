using PhaseDiagram
using Test
using Polymer
using CubicHermiteSpline
using Polyorder
using Scattering

using DelimitedFiles

function readfile(infile)
    data = readdlm(infile, '\t', Float64, '\n'; comments=true)
    # phih  F    mu
    return data[:, 1], data[:, 2], data[:, 3]
end

include("test_surrogate.jl")
include("test_phase.jl")
include("test_diagram.jl")
include("test_phasemodel.jl")
include("test_dis.jl")
include("test_free_energy.jl")
include("test_microphase.jl")
include("test_macrophase.jl")