module Surrogates

using PyCall
using CubicHermiteSpline
using DelaunayTriangulation
using NaturalNeighbours

export scipy_linear_surrogate,
    scipy_cubic_surrogate,
    scipy_nearest_surrogate,
    scipy_kriging_surrogate,
    scipy_chs_surrogate,
    scipy_idw_surrogate,
    scipy_natural_neighbor_surrogate,
    julia_chs_1d_surrogate,
    julia_bchs_2d_surrogate,
    julia_natural_neighbor_surrogate

# 定义全局 Python 库引用
const si = PyNULL()
const pk = PyNULL()
const np = PyNULL()
const sp_dist = PyNULL()
const metpy_interp = PyNULL()

function __init__()
    # 在模块加载时初始化 Python 库
    copy!(si, pyimport("scipy.interpolate"))
    copy!(np, pyimport("numpy"))
    copy!(sp_dist, pyimport("scipy.spatial.distance"))
    try
        copy!(metpy_interp, pyimport("metpy.interpolate"))
    catch
        @warn "metpy 未安装，Natural Neighbor 插值将不可用。请在 Julia 环境中执行: Conda.add(\"metpy\", channel=\"conda-forge\")"
    end
    try
        copy!(pk, pyimport("pykrige.ok"))
    catch
        @warn "pykrige 未安装，Kriging 插值将不可用。请在 Julia 环境中执行: Conda.add(\"pykrige\")"
    end
end

# 辅助函数：通过中心差分计算梯度的闭包
function build_fd_gradients(f; h=1e-5)
    grad_x = (qx, qy) -> (f(qx + h, qy) - f(qx - h, qy)) / (2.0 * h)
    grad_y = (qx, qy) -> (f(qx, qy + h) - f(qx, qy - h)) / (2.0 * h)
    return grad_x, grad_y
end

# ==========================================
# 1. Scipy 封装 (基于 PyCall)
# ==========================================

"""
    scipy_linear_surrogate(coords, values; grad=false)
基于 scipy.interpolate.LinearNDInterpolator 封装的线性插值。
"""
function scipy_linear_surrogate(coords::AbstractMatrix, values::AbstractVector; grad::Bool=false)
    interp_obj = si.LinearNDInterpolator(coords, values)
    f = (x, y) -> interp_obj(x, y)[1]
    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

"""
    scipy_cubic_surrogate(coords, values; grad=false)
基于 scipy.interpolate.CloughTocher2DInterpolator 的三次插值。
"""
function scipy_cubic_surrogate(coords::AbstractMatrix, values::AbstractVector; grad::Bool=false)
    interp_obj = si.CloughTocher2DInterpolator(coords, values)
    f = (x, y) -> interp_obj(x, y)[1]
    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

"""
    scipy_chs_surrogate(coords, values; grad=false)
Scipy 版本的 2D CHS 类似实现 (CloughTocher)。
"""
function scipy_chs_surrogate(coords::AbstractMatrix, values::AbstractVector; grad::Bool=false)
    return scipy_cubic_surrogate(coords, values; grad=grad)
end

"""
    scipy_nearest_surrogate(coords, values; grad=false)
Scipy 最近邻插值。
"""
function scipy_nearest_surrogate(coords::AbstractMatrix, values::AbstractVector; grad::Bool=false)
    interp_obj = si.NearestNDInterpolator(coords, values)
    f = (x, y) -> interp_obj(x, y)[1]
    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

"""
    scipy_kriging_surrogate(coords, values; grad=false)
基于 pykrige 的普通克里金 (Ordinary Kriging) 插值。
"""
function scipy_kriging_surrogate(coords::AbstractMatrix, values::AbstractVector; grad::Bool=false)
    if pk == PyNULL()
        error("pykrige 未能正确加载，无法使用 Kriging 代理模型。")
    end

    x = coords[:, 1]
    y = coords[:, 2]

    ok = pk.OrdinaryKriging(
        x, y, values,
        variogram_model="gaussian",
        variogram_parameters=Dict("sill" => np.var(values), "range" => 1.0, "nugget" => 0.0),
        exact_values=true
    )

    f = (qx, qy) -> begin
        z, ss = ok.execute("points", [qx], [qy])
        return z[1]
    end

    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

"""
    scipy_idw_surrogate(coords, values; p=2, grad=false)
反距离加权插值 (IDW, 距离幂次 p=2)，基于 python numpy 和 cdist 计算。
"""
function scipy_idw_surrogate(coords::AbstractMatrix, values::AbstractVector; p=2, grad::Bool=false)
    x_train = coords[:, 1]
    y_train = coords[:, 2]

    # Return a closure function that acts as a surrogate
    f = (qx, qy) -> begin
        # For single point evaluation [1x2] matrix
        q_coords = [qx qy]
        train_coords = [x_train y_train]

        dist = sp_dist.cdist(q_coords, train_coords)

        # Prevent division by zero
        dist = np.where(dist == 0, 1e-10, dist)

        # Calculate weights 
        weights = 1.0 ./ (dist .^ p)

        # Weighted average Calculation
        z_weights_sum = np.sum(weights .* values', axis=1)
        weight_sum = np.sum(weights, axis=1)

        return (z_weights_sum./weight_sum)[1]
    end

    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

"""
    scipy_natural_neighbor_surrogate(coords, values; grad=false)
自然邻域插值 (Natural Neighbor) 基于 metpy.interpolate.natural_neighbor_to_grid。
"""
function scipy_natural_neighbor_surrogate(coords::AbstractMatrix, values::AbstractVector; grad::Bool=false)
    if metpy_interp == PyNULL()
        error("metpy 未能正确加载，无法使用 Natural Neighbor 代理模型。")
    end

    x_train = coords[:, 1]
    y_train = coords[:, 2]

    f = (qx, qy) -> begin
        # 检查是否直接命中了训练样本点（为了保证严格插值精确穿过样本点）
        idx = findfirst(i -> isapprox(x_train[i], qx, atol=1e-12) && isapprox(y_train[i], qy, atol=1e-12), 1:length(x_train))
        if idx !== nothing
            return values[idx]
        end

        # natural_neighbor_to_grid requires evaluating on grids.
        # Add a tiny perturbation to avoid ZeroDivisionError when grid point matches sample point perfectly. Since we already handled exactly matching above, we do not need the +1e-10 anymore, but we keep it small to protect edges.
        _qx = np.array([qx + 1e-10])
        _qy = np.array([qy + 1e-10])

        # Calling python function
        res = metpy_interp.natural_neighbor_to_grid(x_train, y_train, values, _qx, _qy)

        # In this mode metpy directly returns the interpolated data structure.
        # If it returns a tuple, grab the 3rd element; if just the array, use it directly.
        if isa(res, Tuple)
            return res[3][1]
        else
            return res[1]
        end
    end

    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

# ==========================================
# 2. Julia 原生封装 (基于 CubicHermiteSpline.jl)
# ==========================================

"""
    julia_chs_1d_surrogate(x, y, dydx; grad=false)
基于 CubicHermiteSpline.jl 的一维三次 Hermite 插值。
"""
function julia_chs_1d_surrogate(x::AbstractVector, y::AbstractVector, dydx::AbstractVector; grad::Bool=false)
    spl = UnivariateCHSInterpolation(x, y, dydx)
    f = v -> spl(v)
    if grad
        df = v -> spl(v, grad=true)
        return df
    end
    return f
end

"""
    julia_bchs_2d_surrogate(coords, z, dzdx, dzdy; grad=false)
基于 CubicHermiteSpline.jl 的二维 (Bivariate) 三次 Hermite 插值。
coords: N x 2 矩阵
z: 函数值
dzdx: 对 x 的偏导
dzdy: 对 y 的偏导
"""
function julia_bchs_2d_surrogate(coords::AbstractMatrix, z::AbstractVector, dzdx::AbstractVector, dzdy::AbstractVector; grad::Bool=false)
    x = coords[:, 1]
    y = coords[:, 2]
    spl = BivariateCHSInterpolation(x, y, z, dzdx, dzdy)

    f = (qx, qy) -> spl(qx, qy)

    if grad
        grad_x = (qx, qy) -> begin
            grad_tuple = CubicHermiteSpline.grad(spl, qx, qy)
            return grad_tuple[1]
        end
        grad_y = (qx, qy) -> begin
            grad_tuple = CubicHermiteSpline.grad(spl, qx, qy)
            return grad_tuple[2]
        end
        return grad_x, grad_y
    end

    return f
end

"""
    julia_natural_neighbor_surrogate(coords, z; grad=false)
基于 NaturalNeighbours.jl 的二维自然邻域插值 (Natural Neighbor) 原生实现。
coords: N x 2 矩阵
z: 函数值
"""
function julia_natural_neighbor_surrogate(coords::AbstractMatrix, z::AbstractVector; grad::Bool=false)
    x = coords[:, 1]
    y = coords[:, 2]
    pts = hcat(x, y)' # 2×N matrix expected by some DelaunayTriangulation routines, or simply points

    # 构造 Delaunay 三角剖分
    tri = triangulate(pts)

    # 构造 NaturalNeighbours 插值器
    itp = interpolate(tri, z)

    f = (qx, qy) -> begin
        # 内部求值，NaturalNeighbours 的 interpolate 函数调用可以直接给定 (x, y)
        return itp(qx, qy)
    end

    if grad
        grad_x, grad_y = build_fd_gradients(f)
        return grad_x, grad_y
    end
    return f
end

end # module
