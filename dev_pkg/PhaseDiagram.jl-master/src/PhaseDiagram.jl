module PhaseDiagram

using Reexport
using Setfield

using CubicHermiteSpline
using Polymer
using Scattering
@reexport using Polyorder

using ArgCheck
using LaTeXStrings
using Optim
using NLPModels, NLPModelsModifiers, JSOSolvers
using Roots
using DataInterpolations

include("dis.jl")
import .DIS

include("parameters.jl")
export
    AbstractConfig,
    PolyorderConfig,
    PolymerConfig
export
    getparam,
    setparam!

include("types.jl")
export
    AbstractSurrogate,
    CHSplineSurrogate,
    SplineSurrogate
export
    AbstractEnsemble,
    CanonicalEnsemble,
    GrandCanonicalEnsemble,
    GibbsEnsemble,
    MicrocanonicalEnsemble
export
    AbstractPhase,
    Phase,
    AbstractPhaseDiagram,
    BinaryPhaseDiagram1D,
    BinaryPhaseDiagram2D,
    BinaryPhaseDiagram3D,
    TernaryPhaseDiagram
export
    PolyorderModel,
    polymer_system
export
    AbstractPhaseModel,
    best_model,
    cache_model!,
    MacrophaseModel,
    SurrogatePhaseModel,
    SurrogateMacrophaseModel
export
    MicrophaseModel,
    SurrogateMicrophaseModel
export
    PhaseSeparationType,
    MicrophaseSeparation,
    MacrophaseSeparation
export
    phase_separation_type
export
    DISPhase,
    LAMPhase,
    HEXPhase,
    GYRPhase,
    BCCPhase,
    FCCPhase,
    SCPhase,
    O70Phase,
    SigmaPhase,
    A15Phase,
    DEFAULT_PHASES
export
    description,
    as_variable_name,
    as_ascii_label,
    as_plot_label

include("free_energy.jl")
export
    free_energy,
    gibbs_free_energy,
    grand_free_energy

include("optimize.jl")
export
    OptimOptimizer,
    JSOOptimizer,
    RootsOptimizer,
    OPSOptimizer

include("microphase.jl")
export microphase

include("macrophase.jl")
export macrophase

# include("precompile.jl")

end # module
