abstract type AbstractPhaseDiagram end

abstract type BinaryPhaseDiagram <: AbstractPhaseDiagram end

struct BinaryPhaseDiagram1D{M, X, S<:AbstractParameter{M,X}, P<:AbstractPhase} <: BinaryPhaseDiagram
    xaxis::S
    x::AbstractVector{X}
    phases::AbstractVector{P}
end

struct BinaryPhaseDiagram2D{M, N, X, Y, S<:AbstractParameter{M,X}, T<:AbstractParameter{N,Y}, P<:AbstractPhase} <: BinaryPhaseDiagram
    xaxis::S
    yaxis::T
    x::AbstractVector{X}
    y::AbstractVector{Y}
    phases::AbstractVector{P}
end

struct BinaryPhaseDiagram3D{M, N, O, X, Y, Z, S<:AbstractParameter{M,X}, T<:AbstractParameter{N,Y}, U<:AbstractParameter{O,Z}, P<:AbstractPhase} <: BinaryPhaseDiagram
    xaxis::S
    yaxis::T
    zaxis::U
    x::AbstractVector{X}
    y::AbstractVector{Y}
    z::AbstractVector{Z}
    phases::AbstractVector{P}
end

struct TernaryPhaseDiagram{M, N, O, X, Y, Z, S<:AbstractParameter{M,X}, T<:AbstractParameter{N,Y}, U<:AbstractParameter{O,Z}, P<:AbstractPhase} <: AbstractPhaseDiagram
    leftaxis::S
    rightaxis::T
    bottomaxis::U
    left::AbstractVector{X}
    right::AbstractVector{Y}
    bottom::AbstractVector{Z}
    phases::AbstractVector{P}
end