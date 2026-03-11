abstract type AbstractSurrogate end

"""
    CHSplineSurrogate <: AbstractSurrogate

The input data is interpolated by the Cubic Hermite Spline Interpolation. The corresponding gradients at each data points are required.

This surrogate model is especially useful for multi-component (>=2) polymer system because in such cases the graident is known accompanied by canonical SCFT calculations.
"""
struct CHSplineSurrogate <: AbstractSurrogate
    x
    y
    spl::CubicHermiteSplineInterpolation
end

function CHSplineSurrogate(xs::Vector, ys::Vector, grads::AbstractVector)
    return CHSplineSurrogate(collect(xs), collect(ys), CubicHermiteSplineInterpolation(xs, ys, grads))
end

(chs::CHSplineSurrogate)(x::Number; grad=false) = chs.spl(x; grad=grad)

struct SplineSurrogate{T<:DataInterpolations.AbstractInterpolation} <: AbstractSurrogate
    x
    y
    spl::T
end

function SplineSurrogate(xs, ys; k=3)
    xs = collect(xs)
    ys = collect(ys)
    if k == 1
        return SplineSurrogate(xs, ys, DataInterpolations.LinearInterpolation(ys, xs))
    end
    if k== 2
        return SplineSurrogate(xs, ys, DataInterpolations.QuadraticSpline(ys, xs))
    end
    if k==3
        return SplineSurrogate(xs, ys, DataInterpolations.CubicSpline(ys, xs))
    end
end

(ss::SplineSurrogate)(x::Number) = ss.spl(x)