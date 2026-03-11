# PhaseDiagram.jl [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://liuyxpp.github.io/PhaseDiagram.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://liuyxpp.github.io/PhaseDiagram.jl/dev) [![Build Status](https://github.com/liuyxpp/PhaseDiagram.jl/workflows/CI/badge.svg)](https://github.com/liuyxpp/PhaseDiagram.jl/actions) [![Build Status](https://travis-ci.com/liuyxpp/PhaseDiagram.jl.svg?branch=master)](https://travis-ci.com/liuyxpp/PhaseDiagram.jl) [![Build Status](https://ci.appveyor.com/api/projects/status/github/liuyxpp/PhaseDiagram.jl?svg=true)](https://ci.appveyor.com/project/liuyxpp/PhaseDiagram-jl) [![Coverage](https://codecov.io/gh/liuyxpp/PhaseDiagram.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/liuyxpp/PhaseDiagram.jl)

**PhaseDiagram.jl** is a Julia library for computing phase diagrams of block copolymers and their blends and solutions. For two- or multi-component systems, coexistence boundaries in addition to the phase boundaries will be computed. This package will ultimately replace `PyDiagram`.

## Features

* Computing phase diagrams for macrophase separations of polymer blends and solutions.
* Computing phase diagrams for microphase separations of polymers.

## Usage

### Wrap the SCFT solver

The SCFT solver should be first wrapped. At present, we only have one wrapper, `PolyorderModel`, for wrapping the [Polyorder](https://github.com/liuyxpp/Polyorder.jl) solver. To create a wrapper, we do

```julia
julia> using Polymer  # for AB_A_system
julia> using Scattering  # for UnitCell and BravaisLattice
julia> using Polyorder  # for NoncyclicChainSCFT
julia> using PhaseDiagram  # for PolyorderModel
julia> system = AB_A_system(; χN=16.0);
julia> uc = UnitCell(1.0);
julia> lattice = BravaisLattice(uc);
julia> scft = NoncyclicChainSCFT(system, lattice, 0.1);
# Construct SCFT model, using Polyorder
julia> model = PolyorderModel(scft);
```

### Macrophase separation

To compute the coexistent (binoal) points, we first construct a `MacrophaseModel` and using `macrophase` to solve it.

```julia
# Construct a MacrophaseModel.
# By default, DIS-DIS coexisting point is computed.
julia> mpmodel = MacrophaseModel(model, model; ϕ₁=0.1, ϕ₂=0.7, ϕ₀=0.4);
# Compute the coexistent (binodal) points for this specific `system`.
julia> ϕ₁, ϕ₂, _ = macrophase(OptimOptimizer(), mpmodel);
# Confirm the results are correct.
julia> @test ϕ₁ ≈ 0.05957 atol=1e-4
julia> @test ϕ₂ ≈ 0.84906 atol=1e-4
```

Please see more examples in `test/test_macrophase.jl`.

### Microphase separation

Similarly, to compute the microphase separation, we first construct a `MicrophaseModel` and using `microphase` to solve it.

```julia
julia> fc = fControlParameter(1, 1, fParam)
julia> mipm = MicrophaseModel(model_lam, model_hex, fc, a=0.36, b=0.48, phase1=LAMPhase, phase2=HEXPhase, cached=true)
julia> algo = OPSOptimizer([0.36, 0.42, 0.48], [0.36, 0.42, 0.48])
julia> f, trace = microphase(ops, mipm)
```

Please see more examples in `test/test_microphase.jl`.

## Planned Features

* Adaptive algorithm for constructing phase diagrams.
* Automatic construction of phase diagrams.
* Interactive plots of phase diagrams.
* Publishable plots of phase diagrams.

## Contribute

* Star the package on [github.com](https://github.com/liuyxpp/PhaseDiagram.jl).
* File an issue or make a pull request on [github.com](https://github.com/liuyxpp/PhaseDiagram.jl).
* Contact the author via email <lyx@fudan.edu.cn>.

## Links

* Source code hosted at [github.com](https://github.com/liuyxpp/PhaseDiagram.jl)