module Predictors

using Roots

export predict_phase_equilibrium

function predict_phase_equilibrium(surr_F, surr_mu, x_range)
    # 在连续空间扫描寻找 surrogate 的 spinodal 点
    mus = [surr_mu(p) for p in x_range]
    diff_mus = diff(mus)

    sp1_idx, sp2_idx = -1, -1
    for i in 1:length(diff_mus)-1
        if diff_mus[i] > 0 && diff_mus[i+1] <= 0
            sp1_idx = i + 1
            break
        end
    end
    for i in length(diff_mus)-1:-1:1
        if diff_mus[i] < 0 && diff_mus[i+1] >= 0
            sp2_idx = i + 1
            break
        end
    end

    pa_pred, pb_pred, mu_pred, Fg_pred = NaN, NaN, NaN, NaN

    if sp1_idx != -1 && sp2_idx != -1 && sp1_idx < sp2_idx
        ϕ_sp1_surr = x_range[sp1_idx]
        ϕ_sp2_surr = x_range[sp2_idx]
        m_max = mus[sp1_idx]
        m_min = mus[sp2_idx]

        m_lower_bound = max(m_min, surr_mu(x_range[1]))
        m_upper_bound = min(m_max, surr_mu(x_range[end]))

        if m_upper_bound > m_lower_bound
            function ϕ_alpha_surr(m)
                find_zero(p -> surr_mu(p) - m, (x_range[1], ϕ_sp1_surr), Bisection())
            end
            function ϕ_beta_surr(m)
                find_zero(p -> surr_mu(p) - m, (ϕ_sp2_surr, x_range[end]), Bisection())
            end

            function dF_surr(m)
                pa = ϕ_alpha_surr(m)
                pb = ϕ_beta_surr(m)
                (surr_F(pa) - m * pa) - (surr_F(pb) - m * pb)
            end

            m_tests = range(m_lower_bound + 1e-6, m_upper_bound - 1e-6, length=100)
            dF_vals = Float64[]
            for mt in m_tests
                try
                    push!(dF_vals, dF_surr(mt))
                catch
                    push!(dF_vals, NaN)
                end
            end

            bracket = nothing
            for i in 1:length(dF_vals)-1
                if !isnan(dF_vals[i]) && !isnan(dF_vals[i+1]) && dF_vals[i] * dF_vals[i+1] <= 0
                    bracket = (m_tests[i], m_tests[i+1])
                    break
                end
            end

            if bracket !== nothing
                try
                    mu_pred = find_zero(dF_surr, bracket, Bisection())
                    pa_pred = ϕ_alpha_surr(mu_pred)
                    pb_pred = ϕ_beta_surr(mu_pred)
                    Fg_pred = surr_F(pa_pred) - mu_pred * pa_pred
                catch
                end
            end
        end
    end

    return (ϕ_α=pa_pred, ϕ_β=pb_pred, μ=mu_pred, Fg=Fg_pred)
end

end
