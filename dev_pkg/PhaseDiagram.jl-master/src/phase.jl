"""
There is a great implementation of space group and crystal symmetry in the C library [spglib](https://github.com/spglib/spglib).
It also provide a systematic way to determine the space group given the positions of each atom (or structural unit as in ordered phases in block copolymers).
A binding to Julia is available [spglib.jl](https://github.com/mdavezac/spglib.jl), but it imports partial functions from spglib.
We shall either improve the spglib.jl or in the long-term re-implement it in pure Julia.
Should serve a great project for Student learning crystallography, scattering, and Julia programming.

2022.5.12: Space groups are now supported in Scattering.jl by using Crystalline.jl.
"""

abstract type AbstractPhase{T} end

"""
## Fields
`symmetry`: currently not implemented. Will be implemented in Scattering.jl.
"""
struct Phase{T} <: AbstractPhase{T}
    description::String
    variable_name::String
    ascii_label::String
    plot_label::LaTeXString
    symmetry::Union{SpaceGroup, Nothing}
    function Phase{T}(desc; symmetry=nothing) where T
        varname = String(T)
        ascii_label = ""
        for c in varname
            ascii_label *= unicodesymbol2string(c)
        end
        plot_label = varname
        new{T}(desc, varname, ascii_label, plot_label, symmetry)
    end
end

const DISPhase = Phase{:DIS}("Disordered phase or homogeneous phase")
const LAMPhase = Phase{:LAM}("Lamellar phase"; symmetry=spacegroup(1,1))
const HEXPhase = Phase{:HEX}("Hexagonal cylinder phase"; symmetry=spacegroup(17,2))
const GYRPhase = Phase{:GYR}("Bicontinuous Gyroid phase"; symmetry=spacegroup(230,3))
const BCCPhase = Phase{:BCC}("BCC phase"; symmetry=spacegroup(229,3))
const FCCPhase = Phase{:FCC}("FCC phase"; symmetry=spacegroup(225,3))
const SCPhase = Phase{:SC}("Simple cubic phase"; symmetry=spacegroup(221,3))
const O70Phase = Phase{:O70}("O70 phase"; symmetry=spacegroup(70, 3))
const SigmaPhase = Phase{:Sigma}("Sigma phase"; symmetry=spacegroup(136,3))
const A15Phase = Phase{:A15}("A15 phase"; symmetry=spacegroup(223,3))

const DEFAULT_PHASES = Dict(
    :DIS => DISPhase,
    :LAM => LAMPhase,
    :HEX => HEXPhase,
    :GYR => GYRPhase,
    :BCC => BCCPhase,
    :FCC => FCCPhase,
    :SC => SCPhase,
    :O70 => O70Phase,
    :Sigma => SigmaPhase,
    :A15 => A15Phase,
)

# Shorthand for phase type
const DISType = DISPhase |> typeof
const LAMType = LAMPhase |> typeof
const HEXType = HEXPhase |> typeof
const GYRType = GYRPhase |> typeof
const BCCType = BCCPhase |> typeof
const FCCType = FCCPhase |> typeof
const SCType = SCPhase |> typeof
const O70Type = O70Phase |> typeof
const SigmaType = SigmaPhase |> typeof
const A15Type = A15Phase |> typeof

phase(s::Symbol) = get(DEFAULT_PHASES, s, DISPhase)
phase(s::AbstractString) = phase(Symbol(s))
phasetype(p::T) where {T<:AbstractPhase} = T
phasetype(s) = phasetype(phase(s))

Polymer.description(p::AbstractPhase) = p.description
Polymer.as_variable_name(p::AbstractPhase) = p.variable_name
Polymer.as_ascii_label(p::AbstractPhase) = p.ascii_label
Polymer.as_plot_label(p::AbstractPhase) = p.plot_label