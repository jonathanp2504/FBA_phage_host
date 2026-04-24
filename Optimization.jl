using Statistics

# ============================================================
#  Optimizer voor Benzonase-productie (main branch)
#  Methode: Grid search + Nelder-Mead (pure Julia)
#  Beginbiomassa vastgezet op 1e9
#  Optimaliseert: MOI + infectietijdstip
# ============================================================

const FIXED_BIOMASS = 1e9

function evaluate_benzonase(moi::Float64, t_inf::Float64, p_base::Parameters)
    p_test = Parameters(
        p_base.duration,
        FIXED_BIOMASS,
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
        moi * FIXED_BIOMASS,
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
        @warn "Simulatie gefaald voor MOI=$moi, t_inf=$t_inf: $e"
        return 0.0
    end
end

function run_optimization(p_initial::Parameters, tspan)

    # ----------------------------------------------------------
    # STAP 1: GRID SEARCH
    # ----------------------------------------------------------
    moi_values   = [0.0001, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0]
    t_inf_values = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0, 15.0]

    best_val  = -Inf
    best_moi  = moi_values[1]
    best_tinf = t_inf_values[1]
    total     = length(moi_values) * length(t_inf_values)
    counter   = 0

    println("Grid search gestart ($total combinaties)...")
    println("Biomassa vast op: $FIXED_BIOMASS cellen/L")
    for moi in moi_values
        for t_inf in t_inf_values
            counter += 1
            val = evaluate_benzonase(moi, t_inf, p_initial)
            println("  [$counter/$total] MOI=$moi | t_inf=$t_inf | Benz=$val")
            if val > best_val
                best_val  = val
                best_moi  = moi
                best_tinf = t_inf
            end
        end
    end

    println("\nBeste grid-punt: MOI=$best_moi | t_inf=$best_tinf | Benz=$best_val")

    # ----------------------------------------------------------
    # STAP 2: NELDER-MEAD VERFIJNING
    # ----------------------------------------------------------
    println("\nNelder-Mead verfijning gestart...")

    moi_lb,  moi_ub  = 1e-5, 2.0
    tinf_lb, tinf_ub = 0.5, p_initial.duration - 1.0

    clamp_x(x) = [clamp(x[1], moi_lb,  moi_ub),
                   clamp(x[2], tinf_lb, tinf_ub)]

    f(x) = -evaluate_benzonase(x[1], x[2], p_initial)

    x0    = [best_moi, best_tinf]
    delta = [best_moi * 0.3 + 0.001, 1.0]

    simplex = [
        clamp_x(x0),
        clamp_x(x0 .+ [delta[1], 0.0]),
        clamp_x(x0 .+ [0.0, delta[2]])
    ]
    scores = [f(v) for v in simplex]

    alpha_nm = 1.0
    gamma_nm = 2.0
    rho_nm   = 0.5
    sigma_nm = 0.5
    max_iter = 200
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
                simplex[end] = xe; scores[end] = score_e
            else
                simplex[end] = xr; scores[end] = score_r
            end
        elseif score_r < scores[end]
            simplex[end] = xr; scores[end] = score_r
        else
            xc      = clamp_x(cog .+ rho_nm .* (simplex[end] .- cog))
            score_c = f(xc)
            if score_c < scores[end]
                simplex[end] = xc; scores[end] = score_c
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

    best_x          = simplex[1]
    best_moi_final  = best_x[1]
    best_tinf_final = best_x[2]
    best_benz_final = -scores[1]

    println("\n=== Optimalisatie voltooid ===")
    println("Biomassa (vast)    : $FIXED_BIOMASS cellen/L")
    println("Beste MOI          : $(round(best_moi_final,  digits=4))")
    println("Beste infectietijd : $(round(best_tinf_final, digits=2)) uur")
    println("Max Benzonase      : $(round(best_benz_final, digits=6)) mmol/L")

    return (best_moi       = best_moi_final,
            best_t_inf     = best_tinf_final,
            best_biomass   = FIXED_BIOMASS,
            best_benzonase = best_benz_final)
end