using Statistics
using SciMLBase
include("optimizer_utils.jl")
include("Model4.jl")
# ============================================================
#  OPTIMIZER B — MODEL 4
#  Model 4 heeft b_0 en k_b in plaats van b (variabele burst size)
#  Vaste MOI, optimaliseert infectietijdstip + beginbiomassa
# ============================================================

const FIXED_MOI_B    = 2.0
const W_BENZ_DEFAULT = 0.8
const W_P_DEFAULT    = 0.2

BENZ_REF_B = Ref(1.0)

function setup_optimizer_B(benz_model,
                            exchange_ids::Vector{String},
                            V_max::Vector{Float64},
                            biomass_id::String,
                            benz_id::String,
                            rho_DW::Float64,
                            N0_ref::Float64,
                            duration::Float64)
    println("=== Setup Optimizer B (Model 4): berekening BENZ_REF ===")
    BENZ_REF_B[] = bereken_benz_ref_fba(
        benz_model, exchange_ids, V_max,
        biomass_id, benz_id,
        rho_DW, N0_ref, duration)
end

function evaluate_B(t_inf::Float64, biomass::Float64,
                    p_base::Parameters;
                    w_benz::Float64 = W_BENZ_DEFAULT,
                    w_P::Float64    = W_P_DEFAULT)

    infectiedosis = FIXED_MOI_B * biomass

    p_test = Parameters(
        p_base.duration,
        biomass,
        p_base.alpha_syn,
        p_base.beta_deg,
        p_base.K_s,
        p_base.V_max,
        p_base.p_pref,
        p_base.tau,
        p_base.b_0,
        p_base.k_b,
        p_base.E_coli_cellDW,
        p_base.MW,
        p_base.h_release,
        p_base.biomass_id,
        p_base.ex_ids,
        p_base.all_exchanges,
        p_base.essentials,
        p_base.fbaModelNaive,
        p_base.fbaModelLysogen,
        p_base.mu_N,
        copy(p_base.q_N),
        p_base.mu_l,
        copy(p_base.q_l),
        p_base.q_benz_l,
        p_base.k_attach,
        p_base.k_dettach,
        p_base.k_inject,
        p_base.K_mal,
        t_inf,
        infectiedosis,
        p_base.benz_id,
        p_base.k_tox,
        p_base.beta_benz,
        p_base.q_benz,
        copy(p_base.mu_max),
        copy(p_base.e_max),
        p_base.f_prod
    )

    try
        sol = run(p_test)
        !SciMLBase.successful_retcode(sol) && return -Inf

        max_benz  = maximum(sol[Benzind, :])
        benz_norm = max_benz / BENZ_REF_B[]

        P_max_sim = maximum(sol[Pfind, :])
        P_norm    = P_max_sim > 0.0 ? infectiedosis / P_max_sim : 0.0

        return w_benz * benz_norm - w_P * P_norm

    catch e
        @warn "Simulatie gefaald t_inf=$t_inf, biomass=$biomass: $e"
        return -Inf
    end
end

function run_optimization_B(p_initial::Parameters;
                             w_benz::Float64 = W_BENZ_DEFAULT,
                             w_P::Float64    = W_P_DEFAULT)

    t_inf_values   = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]
    biomass_values = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]

    best_score   = -Inf
    best_tinf    = t_inf_values[2]
    best_biomass = biomass_values[3]
    total        = length(t_inf_values) * length(biomass_values)
    counter      = 0

    println("\nGrid search B Model 4 ($total combinaties) | MOI vast = $FIXED_MOI_B | w_benz=$w_benz")
    for t_inf in t_inf_values
        for biomass in biomass_values
            counter += 1
            score = evaluate_B(t_inf, biomass, p_initial; w_benz=w_benz, w_P=w_P)
            @printf("  [%3d/%d] t_inf=%-5.1f | biomass=%-10.2e | score=%.5f\n",
                counter, total, t_inf, biomass, isfinite(score) ? score : -999.0)
            if isfinite(score) && score > best_score
                best_score = score; best_tinf = t_inf; best_biomass = biomass
            end
        end
    end
    println("\nBeste grid-punt: t_inf=$best_tinf | biomass=$best_biomass | score=$(round(best_score,digits=5))")

    println("\nNelder-Mead verfijning B (Model 4)...")
    f_nm(x) = evaluate_B(x[1], x[2], p_initial; w_benz=w_benz, w_P=w_P)
    x0    = [best_tinf, best_biomass]
    delta = [1.0, best_biomass * 0.3]
    lb    = [0.5, 1e4]
    ub    = [p_initial.duration - 1.0, 1e11]

    best_x, best_score_nm = nelder_mead(f_nm, x0, delta, lb, ub)

    println("\n=== Optimalisatie B (Model 4) voltooid ===")
    println("MOI (vast)         : $FIXED_MOI_B")
    @printf("Infectietijdstip   : %.2f h\n",       best_x[1])
    @printf("Beginbiomassa      : %.3e cellen/L\n", best_x[2])
    @printf("Infectiedosis      : %.3e fagen/L\n",  FIXED_MOI_B * best_x[2])
    @printf("Gewogen score      : %.6f\n",           best_score_nm)

    return (best_moi   = FIXED_MOI_B,
            best_t_inf = best_x[1],
            best_N0    = best_x[2],
            best_score = best_score_nm)
end

function run_pareto_B(p_initial::Parameters;
                      n_weights::Int  = 11,
                      figname::String = "pareto_B_M4.png")

    function single_opt(w::Float64)
        res           = run_optimization_B(p_initial; w_benz=w, w_P=1.0-w)
        infectiedosis = FIXED_MOI_B * res.best_N0
        p_eval = Parameters(
            p_initial.duration, res.best_N0,
            p_initial.alpha_syn, p_initial.beta_deg,
            p_initial.K_s, p_initial.V_max, p_initial.p_pref,
            p_initial.tau, p_initial.b_0, p_initial.k_b,
            p_initial.E_coli_cellDW, p_initial.MW, p_initial.h_release,
            p_initial.biomass_id, p_initial.ex_ids,
            p_initial.all_exchanges, p_initial.essentials,
            p_initial.fbaModelNaive, p_initial.fbaModelLysogen,
            p_initial.mu_N, copy(p_initial.q_N),
            p_initial.mu_l, copy(p_initial.q_l), p_initial.q_benz_l,
            p_initial.k_attach, p_initial.k_dettach, p_initial.k_inject,
            p_initial.K_mal,
            res.best_t_inf, infectiedosis,
            p_initial.benz_id, p_initial.k_tox,
            p_initial.beta_benz, p_initial.q_benz,
            copy(p_initial.mu_max), copy(p_initial.e_max), p_initial.f_prod)

        sol      = run(p_eval)
        max_benz = SciMLBase.successful_retcode(sol) ? maximum(sol[Benzind, :]) : 0.0
        P_max    = SciMLBase.successful_retcode(sol) ? maximum(sol[Pfind,   :]) : Inf

        return (w=w, benz=max_benz, P_max=P_max,
                t_inf=res.best_t_inf, moi=FIXED_MOI_B,
                N0=res.best_N0, score=res.best_score)
    end

    return pareto_front(single_opt; n_weights=n_weights, figname=figname)
end
