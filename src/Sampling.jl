module Sampling

export sample_next_points

"""
    sample_next_points(ϕ_pred, current_x; strategy="OPS")

根据 Predicators 输出的预测相平衡点 `ϕ_pred`，计算下一步真正要送入 Black-box 评估的新抽样点集合。

参数:
- `ϕ_pred`: Tuple 或 Vector，即 predictor 给出的平衡点。对于 2组分(1维空间)，通常是 2个值的 tuple (ϕ_α, ϕ_β)。对于多维体系，由于预测点可能是向量，故可以用数组或元组的形式给出多个成分。
- `current_x`: 当前已经被 Black-box 评估过的自变量点集合 (1维时通常是排序好的 array，多维时则是点集数组)。
- `strategy`: 采样策略，支持 "OPS"(On-Position), "BS"(Bisection), "TS"(Tri-section), "QS"(Quad-section), "RS"(Random)。

返回:
- 一个数组，包含需要进行下一轮评估的全新的坐标点。
"""
function sample_next_points(ϕ_pred, current_x; strategy="OPS")
    # 清洗：剔除 NaN 失败的预测，保留有效的预测点
    valid_preds = filter(p -> !any(isnan.(p)), ϕ_pred)
    new_points = similar(valid_preds, 0)
    
    if isempty(valid_preds)
        # 如果 Predictor 返回完全失败，启用全局随机回退机制
        return _fallback_random(current_x)
    end
    
    # 检查输入是 1D 还是 2D 体系（通过预测点的类型，如果是单个标量数字说明在 1 维体系，也就是 2 组分）
    is_1d = eltype(valid_preds) <: Number

    if strategy == "OPS"
        # OPS (On-Position Sampling) : 在预测点直接打点。1D 和 高维通用。
        for vp in valid_preds
            if !_is_too_close(vp, current_x, 1e-7)
                push!(new_points, vp)
            end
        end
        return isempty(new_points) ? _fallback_random(current_x) : new_points
    end
    
    # ================= 以下策略主要针对一维（二组分系统）生效 =================
    if !is_1d
        # 三组分(2D空间)或更高维：在非 OPS 策略时统一降级到 OPS 处理，因为网格分割在高维不适用或太昂贵
        @warn "Strategy '$strategy' is primarily designed for 1D phase boundaries. Falling back to OPS for high-dimensional predictions."
        return sample_next_points(ϕ_pred, current_x, strategy="OPS")
    end

    # 对于 1D 体系特有的细分网格采样逻辑，这需要基于预测点，在其所处的两个相邻采样格子内执行。
    # 因为 valid_preds 可能有多个（如 ϕ_α, ϕ_β 等），我们分别对其操作：
    for p_val in valid_preds
        # 寻找包含 p_val 的最小已知区间 (phi_i, phi_j)
        # 这里假设 current_x 在 1D 时是排序好的标量数组
        idx = searchsortedfirst(current_x, p_val)
        
        # 边界处理：如果预测点落出外推区域，回退为 OPS
        if idx == 1 || idx > length(current_x)
            if !_is_too_close(p_val, current_x, 1e-7)
                push!(new_points, p_val)
            end
            continue
        end

        phi_i = current_x[idx - 1]
        phi_j = current_x[idx]
        interval_len = phi_j - phi_i
        
        if strategy == "BS" # Bisection (中间加1点)
            pt = phi_i + interval_len / 2.0
            if !_is_too_close(pt, current_x, 1e-7)
                push!(new_points, pt)
            end
            
        elseif strategy == "TS" # Tri-section (加2点)
            for k in 1:2
                pt = phi_i + k * interval_len / 3.0
                if !_is_too_close(pt, current_x, 1e-7)
                    push!(new_points, pt)
                end
            end
            
        elseif strategy == "QS" # Quad-section (加3点)
            for k in 1:3
                pt = phi_i + k * interval_len / 4.0
                if !_is_too_close(pt, current_x, 1e-7)
                    push!(new_points, pt)
                end
            end
            
        elseif strategy == "RS" # Random (在区间内随机打一点)
            pt = phi_i + rand() * interval_len
            if !_is_too_close(pt, current_x, 1e-7)
                push!(new_points, pt)
            end
            
        else
            error("Unknown local sampling strategy: $strategy")
        end
    end

    # 如果细分采样提取出的全都是已经被估值的重合点，回退
    if isempty(new_points)
        new_points = sample_next_points(ϕ_pred, current_x, strategy="OPS")
    end
    
    return isempty(new_points) ? _fallback_random(current_x) : new_points
end

# 内部辅导函数
function _is_too_close(pt, X_hist, tol=1e-7)
    if isempty(X_hist)
        return false
    end
    if pt isa Number
        # 1D 系
        return minimum(abs.(X_hist .- pt)) < tol
    else
        # 2D 或 多维
        # 计算欧式距离
        dists = [sqrt(sum((x .- pt).^2)) for x in X_hist]
        return minimum(dists) < tol
    end
end

function _fallback_random(current_x)
    new_points = similar(current_x, 0)
    if isempty(current_x)
        push!(new_points, rand())
        return new_points
    end
    
    if current_x[1] isa Number
        # 1D 均匀产生一个
        push!(new_points, rand() * (maximum(current_x) - minimum(current_x)) + minimum(current_x))
    else
        # 高维：在数据的包围盒中生成一个
        mins = [minimum(v) for v in eachrow(hcat(current_x...))]
        maxs = [maximum(v) for v in eachrow(hcat(current_x...))]
        pt = [rand() * (maxs[i] - mins[i]) + mins[i] for i=1:length(mins)]
        push!(new_points, pt)
    end
    return new_points
end

end # module
