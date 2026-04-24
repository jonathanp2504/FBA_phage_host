using Statistics

# ============================================================
#  Optimizer voor Benzonase-productie
#  Methode: Grid search + Nelder-Mead (pure Julia)
#  MOI vastgezet op 2.0
#  Optimaliseert: infectietijdstip + beginbiomassa
# ============================================================

function evaluate_benzonase(moi::Float64, t_inf::Float64, biomass::Float64, p_base::Parameters)
    p_test = Parameters(
        p_base.duration,
        biomass,
        p_base.alpha_syn,
        p_base.beta_deg,
        p_base.K_s,
        p_base.V_max,
        p_base.p_pref,
        p_base.tau,
        p_base.b,
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
        moi * biomass,
        p_base.benz_id,
        p_base.k_tox,
        p_base.beta_benz,
        p_base.q_benz,
        copy(p_base.mu_max),
        copy(p_base.e_max),
        p_base.f_prod,
        p_base.Y_benz
    )

    try
        sol = run(p_test)
        if SciMLBase.successful_retcode(sol)
            return maximum(sol[Benzind, :])
        else
            return 0.0
        end
    catch e
        @warn "Simulatie gefaald voor MOI=$moi, t_inf=$t_inf, biomass=$biomass: $e"
        return 0.0
    end
end

function run_optimization(p_initial::Parameters, tspan)

    # ----------------------------------------------------------
    # STAP 1: GRID SEARCH
    # MOI is vastgezet op 2.0
    # ----------------------------------------------------------
    moi_values     = [2.0]
    t_inf_values   = [1.0, 2.0, 3.0, 5.0, 7.0]
    biomass_values = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]

    best_val     = -Inf
    best_moi     = 2.0
    best_tinf    = t_inf_values[1]
    best_biomass = biomass_values[1]
    total        = length(moi_values) * length(t_inf_values) * length(biomass_values)
    counter      = 0

    println("Grid search gestart ($total combinaties)...")
    for moi in moi_values
        for t_inf in t_inf_values
            for biomass in biomass_values
                counter += 1
                val = evaluate_benzonase(moi, t_inf, biomass, p_initial)
                println("  [$counter/$total] MOI=$moi | t_inf=$t_inf | biomass=$biomass | Benz=$val")
                if val > best_val
                    best_val     = val
                    best_moi     = moi
                    best_tinf    = t_inf
                    best_biomass = biomass
                end
            end
        end
    end

    println("\nBeste grid-punt: MOI=$best_moi | t_inf=$best_tinf | biomass=$best_biomass | Benz=$best_val")

    # ----------------------------------------------------------
    # STAP 2: NELDER-MEAD VERFIJNING
    # Optimaliseert t_inf en biomassa (MOI blijft 2.0)
    # ----------------------------------------------------------
    println("\nNelder-Mead verfijning gestart...")

    # Grenzen
    tinf_lb,    tinf_ub    = 0.0, p_initial.duration - 1.0
    biomass_lb, biomass_ub = 1e4, 1e10

    # Hulpfunctie: klem waarde binnen grenzen
    clamp_x(x) = [clamp(x[1], tinf_lb, tinf_ub),
                   clamp(x[2], biomass_lb, biomass_ub)]

    # Objectieffunctie (minimaliseren = negatieve Benzonase)
    f(x) = -evaluate_benzonase(2.0, x[1], x[2], p_initial)

    # Initiële simplex rond beste grid-punt
    x0    = [best_tinf, best_biomass]
    delta = [1.0, best_biomass * 0.3]

    simplex = [
        clamp_x(x0),
        clamp_x(x0 .+ [delta[1], 0.0]),
        clamp_x(x0 .+ [0.0, delta[2]])
    ]
    scores = [f(v) for v in simplex]

    # Nelder-Mead parameters
    alpha_nm = 1.0
    gamma_nm = 2.0
    rho_nm   = 0.5
    sigma_nm = 0.5

    max_iter = 100
    tol      = 1e-4

    for iter in 1:max_iter
        ord     = sortperm(scores)
        simplex = simplex[ord]
        scores  = scores[ord]

        if abs(scores[end] - scores[1]) < tol
            println("  Geconvergeerd na $iter iteraties.")
            break
        end

        n   = length(simplex)
        cog = clamp_x(mean(simplex[1:n-1]))

        xr      = clamp_x(cog .+ alpha_nm .* (cog .- simplex[end]))
        score_r = f(xr)

        if score_r < scores[1]
            xe      = clamp_x(cog .+ gamma_nm .* (xr .- cog))
            score_e = f(xe)
            if score_e < score_r
                simplex[end] = xe
                scores[end]  = score_e
            else
                simplex[end] = xr
                scores[end]  = score_r
            end
        elseif score_r < scores[end]
            simplex[end] = xr
            scores[end]  = score_r
        else
            xc      = clamp_x(cog .+ rho_nm .* (simplex[end] .- cog))
            score_c = f(xc)
            if score_c < scores[end]
                simplex[end] = xc
                scores[end]  = score_c
            else
                for i in 2:n
                    simplex[i] = clamp_x(simplex[1] .+ sigma_nm .* (simplex[i] .- simplex[1]))
                    scores[i]  = f(simplex[i])
                end
            end
        end

        if iter % 20 == 0
            println("  Iter $iter | beste Benz=$(round(-scores[1], digits=6))")
        end
    end

    # Eindresultaat
    best_x             = simplex[1]
    best_tinf_final    = best_x[1]
    best_biomass_final = best_x[2]
    best_benz_final    = -scores[1]

    println("\n=== Optimalisatie voltooid ===")
    println("MOI (vast)         : 2.0")
    println("Beste infectietijd : $(round(best_tinf_final, digits=2)) uur")
    println("Beste beginbiomassa: $(round(best_biomass_final, digits=0)) cellen/L")
    println("Max Benzonase      : $(round(best_benz_final, digits=6)) mmol/L")

    return (best_moi       = 2.0,
            best_t_inf     = best_tinf_final,
            best_biomass   = best_biomass_final,
            best_benzonase = best_benz_final)
end