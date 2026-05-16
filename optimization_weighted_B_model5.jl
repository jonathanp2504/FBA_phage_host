using Statistics
using SciMLBase

# ============================================================
#  GEWOGEN OPTIMIZER: 80% Benzonase + 5% minimale P/N
#  Versie B — MODEL 5 SPECIFIEK
#  Model 5 heeft fbaModelLytic, phage_id en q_phage_L
#  als extra velden en heeft geen b veld.
#
#  Optimaliseert infectietijdstip + beginbiomassa
#  MOI vast op FIXED_MOI_B
# ============================================================

const BENZ_REF_B  = 2.0e-3
const PN_REF_B    = 1.06e12
const FIXED_MOI_B = 2.0

function evaluate_weighted_B(t_inf::Float64, biomass::Float64,
                              p_base::Parameters)
    p_test = Parameters(
        p_base.duration,
        biomass,
        p_base.alpha_syn,
        p_base.beta_deg,
        p_base.K_s,
        p_base.V_max,
        p_base.p_pref,
        p_base.tau,
        p_base.E_coli_cellDW,
        p_base.MW,
        p_base.h_release,
        p_base.biomass_id,
        p_base.ex_ids,
        p_base.all_exchanges,
        p_base.essentials,
        p_base.fbaModelNaive,
        p_base.fbaModelLysogen,
        p_base.fbaModelLytic,    # model 5: lytische FBA cache
        p_base.mu_N,
        copy(p_base.q_N),
        p_base.mu_l,
        copy(p_base.q_l),
        p_base.q_benz_l,
        p_base.q_phage_L,        # model 5: faagproductieflux
        p_base.k_attach,
        p_base.k_dettach,
        p_base.k_inject,
        p_base.K_mal,
        t_inf,
        FIXED_MOI_B * biomass,
        p_base.benz_id,
        p_base.phage_id,         # model 5: phage reactie ID
        p_base.k_tox,
        p_base.beta_benz,
        p_base.q_benz,
        copy(p_base.mu_max),
        copy(p_base.e_max),
        p_base.f_prod
    )

    try
        sol = run(p_test)
        if !SciMLBase.successful_retcode(sol)
            return -Inf
        end

        max_benz  = maximum(sol[Benzind, :])
        benz_norm = max_benz / BENZ_REF_B

        final_P  = sol[Pfind, end]
        final_N  = sol[Nind,  end]
        PN_ratio = final_N > 1.0 ? final_P / final_N : 1e12
        PN_norm  = PN_ratio / PN_REF_B

        return 0.8 * benz_norm - 0.2 * PN_norm

    catch e
        @warn "Simulatie gefaald t_inf=$t_inf, biomass=$biomass: $e"
        return -Inf
    end
end

function run_optimization_B(p_initial::Parameters)

    t_inf_values   = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]
    biomass_values = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]

    best_score   = -Inf
    best_tinf    = t_inf_values[1]
    best_biomass = biomass_values[1]
    total        = length(t_inf_values) * length(biomass_values)
    counter      = 0

    println("Grid search B ($total combinaties) | MOI vast = $FIXED_MOI_B")
    for t_inf in t_inf_values
        for biomass in biomass_values
            counter += 1
            score = evaluate_weighted_B(t_inf, biomass, p_initial)
            # Sla scores van exact 0.0 over — geen productie/fagen, biologisch betekenisloos
            if score == 0.0
                println("  [$counter/$total] t_inf=$t_inf | biomass=$biomass | score=0.0 (overgeslagen)")
                continue
            end
            println("  [$counter/$total] t_inf=$t_inf | biomass=$biomass | score=$(round(score, digits=5))")
            if score > best_score
                best_score   = score
                best_tinf    = t_inf
                best_biomass = biomass
            end
        end
    end

    # Als alle grid-punten nul of -Inf zijn
    if best_score == -Inf
        @warn "Alle grid-punten gaven score 0.0 of -Inf. Controleer kalibratie en parameters."
        return (best_moi=FIXED_MOI_B, best_t_inf=t_inf_values[2],
                best_biomass=biomass_values[3], best_score=-Inf)
    end

    println("\nBeste grid-punt: t_inf=$best_tinf | biomass=$best_biomass | score=$best_score")

    println("\nNelder-Mead verfijning B...")
    function f_nm(x)
        score = evaluate_weighted_B(x[1], x[2], p_initial)
        # Behandel score 0.0 als -Inf voor Nelder-Mead
        return score == 0.0 ? -Inf : score
    end

    x0    = [best_tinf, best_biomass]
    delta = [1.0, best_biomass * 0.3]
    lb    = [0.5,  1e4]
    ub    = [p_initial.duration - 1.0, 1e11]

    best_x, best_score_nm = nelder_mead(f_nm, x0, delta, lb, ub)

    println("\n=== Optimalisatie B voltooid ===")
    println("MOI (vast)         : $FIXED_MOI_B")
    println("Infectietijdstip   : $(round(best_x[1],     digits=2)) uur")
    println("Beginbiomassa      : $(round(best_x[2],     digits=0)) cellen/L")
    println("Gewogen score      : $(round(best_score_nm,  digits=6))")

    return (best_moi     = FIXED_MOI_B,
            best_t_inf   = best_x[1],
            best_biomass = best_x[2],
            best_score   = best_score_nm)
end