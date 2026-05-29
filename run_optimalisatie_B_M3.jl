using Statistics
using SciMLBase
using Printf
include("optimizer_utils.jl")
import Main.Model3: run, Benzind, Pfind, Nind

# ============================================================
#  OPTIMIZER B — MODEL 3
#  Vaste MOI (0.01, optimaal uit Optimizer A),
#  optimaliseert infectietijdstip + beginbiomassa
# ============================================================

const FIXED_MOI_B = 0.01

const BENZ_REF_B = Ref(0.124)  # mmol/L, zelfde referentie als Optimizer A

function evaluate_B(t_inf::Float64, biomass::Float64,
                    p_base;
                    w_benz::Float64 = 10.0,
                    w_P::Float64    = 0.1,
                    w_t::Float64    = 0.1)

    infectiedosis = FIXED_MOI_B * biomass

    p_test = Model3.Parameters(
        p_base.duration, biomass,
        p_base.alpha_syn, p_base.beta_deg, p_base.K_s, p_base.V_max, p_base.p_pref,
        p_base.tau, p_base.b, p_base.E_coli_cellDW, p_base.MW, p_base.h_release,
        p_base.biomass_id, p_base.ex_ids, p_base.all_exchanges, p_base.essentials,
        p_base.fbaModelNaive, p_base.fbaModelLysogen,
        p_base.mu_N, copy(p_base.q_N), p_base.mu_l, copy(p_base.q_l), p_base.q_benz_l,
        p_base.k_attach, p_base.k_dettach, p_base.k_inject, p_base.K_mal,
        t_inf, infectiedosis,
        p_base.benz_id, p_base.k_tox, p_base.beta_benz, p_base.q_benz,
        copy(p_base.mu_max), copy(p_base.e_max), p_base.f_prod)

    try
        sol = redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                Model3.run(p_test)
            end
        end
        !SciMLBase.successful_retcode(sol) && return -Inf

        max_benz   = maximum(sol[Benzind, :])
        benz_norm  = 1.0 - exp(-max_benz / BENZ_REF_B[])
        P_max_sim  = maximum(sol[Pfind, :])
        P_norm     = P_max_sim > 0.0 ? 1.0 - exp(-log(infectiedosis) / log(P_max_sim)) : 0.0
        t_max_benz = sol.t[findfirst(sol[Benzind, :] .>= 0.95 * max_benz)]
        t_norm     = 1.0 - exp(-max(0.0, t_max_benz - t_inf) / p_base.duration)

        return w_benz * benz_norm - w_P * P_norm - w_t * t_norm

    catch e
        @warn "Simulatie gefaald t_inf=$t_inf, biomass=$biomass: $e"
        return -Inf
    end
end

function run_optimization_B(p_initial;
                             w_benz::Float64 = 10.0,
                             w_P::Float64    = 0.1,
                             w_t::Float64    = 0.1)

    t_inf_values   = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]
    biomass_values = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]

    best_score   = -Inf
    best_tinf    = t_inf_values[2]
    best_biomass = biomass_values[3]
    total        = length(t_inf_values) * length(biomass_values)
    counter      = 0

    println("\nGrid search B M3 ($total combinaties) | MOI vast = $FIXED_MOI_B | w_benz=$w_benz")
    for t_inf in t_inf_values
        for biomass in biomass_values
            counter += 1
            score = evaluate_B(t_inf, biomass, p_initial; w_benz=w_benz, w_P=w_P, w_t=w_t)
            @printf("  [%3d/%d] t_inf=%-5.1f | biomass=%-10.2e | score=%.5f\n",
                counter, total, t_inf, biomass, isfinite(score) ? score : -999.0)
            if isfinite(score) && score > best_score
                best_score = score; best_tinf = t_inf; best_biomass = biomass
            end
        end
    end
    println("\nBeste grid-punt: t_inf=$best_tinf | biomass=$best_biomass | score=$(round(best_score, digits=5))")

    println("\nNelder-Mead verfijning B (M3)...")
    f_nm(x) = evaluate_B(x[1], x[2], p_initial; w_benz=w_benz, w_P=w_P, w_t=w_t)
    x0    = [best_tinf, best_biomass]
    delta = [1.0, best_biomass * 0.3]
    lb    = [0.5, 1e4]
    ub    = [p_initial.duration - 1.0, 1e11]

    best_x, best_score_nm = nelder_mead(f_nm, x0, delta, lb, ub)

    println("\n=== Optimalisatie B (M3) voltooid ===")
    println("MOI (vast)         : $FIXED_MOI_B")
    @printf("Infectietijdstip   : %.2f h\n",       best_x[1])
    @printf("Beginbiomassa      : %.3e cellen/L\n", best_x[2])
    @printf("Infectiedosis      : %.3e fagen/L\n",  FIXED_MOI_B * best_x[2])
    @printf("Gewogen score      : %.6f\n",           best_score_nm)

    return (best_moi=FIXED_MOI_B, best_t_inf=best_x[1],
            best_N0=best_x[2], best_score=best_score_nm)
end
