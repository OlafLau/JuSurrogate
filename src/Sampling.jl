module Sampling

export sample_next_points

function sample_next_points(pa_pred, pb_pred, current_x; strategy="OPS")
    new_points = Float64[]

    if strategy == "OPS"
        if !isnan(pa_pred) && !isnan(pb_pred)
            if minimum(abs.(current_x .- pa_pred)) > 1e-7
                push!(new_points, pa_pred)
            end
            if minimum(abs.(current_x .- pb_pred)) > 1e-7
                push!(new_points, pb_pred)
            end
        end
        # 如果预测失败，或者预测点与现有数据重合，采用 RS 随机采样后备方案
        if isempty(new_points)
            rand_ϕ = rand() * (current_x[end] - current_x[1]) + current_x[1]
            push!(new_points, rand_ϕ)
        end
    elseif strategy == "RS"
        if !isnan(pa_pred) && !isnan(pb_pred)
            # 在预测点区间内随机采样一个点
            pmin = min(pa_pred, pb_pred)
            pmax = max(pa_pred, pb_pred)
            rand_ϕ = rand() * (pmax - pmin) + pmin
            push!(new_points, rand_ϕ)
        else
            rand_ϕ = rand() * (current_x[end] - current_x[1]) + current_x[1]
            push!(new_points, rand_ϕ)
        end
    end

    return new_points
end

end
