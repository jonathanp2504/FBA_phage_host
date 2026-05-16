using Statistics
using SciMLBase

# ============================================================
#  GEWOGEN OPTIMIZER: 80% Benzonase + 20% minimale P/N
#  Versie A: optimaliseert MOI + infectietijdstip
#            beginbiomassa vast op FIXED_BIOMASS
#
#  Gebruik voor: Model 2, 3, 4, 5
#
#  HOE REFERENTIEWAARDEN INSTELLEN:
#  1. Draai eerst run(p) met je standaardparameters
#  2. Kijk naar maximum(sol[Benzind,:]) -> dit is BENZ_REF
#  3. Kijk naar sol[Pfind,end]/sol[Nind,end] -> dit is PN_REF
#  4. Pas de constanten hieronder aan
# ============================================================

# Referentiewaarden voor normalisatie
# Pas deze aan na een eerste exploratieve run
const BENZ_REF_A  = 2.0e-3   # typische max Benzonase [mmol/L]
const PN_REF_A    = 3.0e9     # typische P/N eindverhouding [-]
const FIXED_BIOMASSA = 1e9  # vaste beginbiomassa [cellen/L]
const FIXED_MOI_A    = 2.0  # vaste MOI

# ============================================================
#  Objectieffunctie
# ============================================================
function evaluate_weighted_A(t_inf::Float64, moi::Float64,
                              biomass::Float64, p_base::Parameters)
    p_test = Parameters(
        p_base.duration,
        biomass,
        p_base.alpha_syn,
        p_base.beta_deg,
        p_base.K_s,
        p_base.V_max,
        p_base.p_pref,
        p_base.tau,
        p_base.b,      #bij model 2,3 b toevoegen, hier eruit
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
        moi * biomass,       # infectiedosis = MOI × N0
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
        if !SciMLBase.successful_retcode(sol)
            return -Inf
        end

        # Term 1: maximale Benzonase productie
        max_benz  = maximum(sol[Benzind, :])
        benz_norm = max_benz / BENZ_REF_A

        # Term 2: minimale P/N verhouding op het einde
        # P/N = vrije fagen / naïeve cellen
        # Lage P/N = weinig vrije fagen in het medium = veiliger voor labo
        final_P = sol[Pfind, end]
        final_N = sol[Nind, end]
        PN_ratio = final_N > 1.0 ? final_P / final_N : 1e12
        PN_norm  = PN_ratio / PN_REF_A

        # Gewogen objectief: maximaliseer (hoog = goed)
        # + benz_norm: hogere Benzonase is beter
        # - PN_norm:   lagere P/N is beter (vandaar het minteken)
        return 0.8 * benz_norm - 0.2 * PN_norm

    catch e
        @warn "Simulatie gefaald t_inf=$t_inf, moi=$moi: $e"
        return -Inf
    end
end

# ============================================================
#  Nelder-Mead kern (herbruikbare functie)
# ============================================================
function nelder_mead(f, x0, delta, lb, ub;
                     max_iter=200, tol=1e-4,
                     alpha=1.0, gamma=2.0, rho=0.5, sigma=0.5)

    clamp_x(x) = [clamp(x[i], lb[i], ub[i]) for i in eachindex(x)]
    n = length(x0)

    # Initialiseer simplex: startpunt + n verstoringen
    simplex = Vector{Vector{Float64}}(undef, n+1)
    simplex[1] = clamp_x(x0)
    for i in 1:n
        v = copy(x0)
        v[i] += delta[i]
        simplex[i+1] = clamp_x(v)
    end
    scores = [f(v) for v in simplex]

    for iter in 1:max_iter
        ord     = sortperm(scores, rev=true)   # hoogste score eerst
        simplex = simplex[ord]
        scores  = scores[ord]

        # Convergentiecheck
        if abs(scores[1] - scores[end]) < tol
            println("  Nelder-Mead geconvergeerd na $iter iteraties.")
            break
        end

        # Zwaartepunt van de beste n punten
        cog = clamp_x(mean(simplex[1:n]))

        # Reflectie
        xr      = clamp_x(cog .+ alpha .* (cog .- simplex[end]))
        score_r = f(xr)

        if score_r > scores[1]
            # Expansie
            xe      = clamp_x(cog .+ gamma .* (xr .- cog))
            score_e = f(xe)
            if score_e > score_r
                simplex[end] = xe; scores[end] = score_e
            else
                simplex[end] = xr; scores[end] = score_r
            end
        elseif score_r > scores[end]
            simplex[end] = xr; scores[end] = score_r
        else
            # Contractie
            xc      = clamp_x(cog .+ rho .* (simplex[end] .- cog))
            score_c = f(xc)
            if score_c > scores[end]
                simplex[end] = xc; scores[end] = score_c
            else
                # Krimp
                for i in 2:n+1
                    simplex[i] = clamp_x(simplex[1] .+ sigma .* (simplex[i] .- simplex[1]))
                    scores[i]  = f(simplex[i])
                end
            end
        end

        if iter % 25 == 0
            println("  Iter $iter | beste score = $(round(scores[1], digits=6))")
        end
    end

    return simplex[1], scores[1]
end

# ============================================================
#  run_optimization_A: optimaliseert t_inf + MOI
#  beginbiomassa vast op FIXED_BIOMASSA
# ============================================================
function run_optimization_A(p_initial::Parameters)

    # --- Grid search ---
    moi_values   = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    t_inf_values = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]

    best_score = -Inf
    best_moi   = moi_values[1]
    best_tinf  = t_inf_values[1]
    total      = length(moi_values) * length(t_inf_values)
    counter    = 0

    println("Grid search A ($total combinaties) | biomassa vast = $FIXED_BIOMASSA")
    for moi in moi_values
        for t_inf in t_inf_values
            counter += 1
            score = evaluate_weighted_A(t_inf, moi, FIXED_BIOMASSA, p_initial)
            println("  [$counter/$total] MOI=$moi | t_inf=$t_inf | score=$(round(score, digits=5))")
            if score > best_score
                best_score = score
                best_moi   = moi
                best_tinf  = t_inf
            end
        end
    end
    println("\nBeste grid-punt: MOI=$best_moi | t_inf=$best_tinf | score=$best_score")

    # --- Nelder-Mead verfijning ---
    println("\nNelder-Mead verfijning A...")
    f_nm(x) = evaluate_weighted_A(x[1], x[2], FIXED_BIOMASSA, p_initial)

    x0    = [best_tinf, best_moi]
    delta = [1.0, best_moi * 0.3 + 1e-4]
    lb    = [0.5,  1e-5]
    ub    = [p_initial.duration - 1.0, 10.0]

    best_x, best_score_nm = nelder_mead(f_nm, x0, delta, lb, ub)

    println("\n=== Optimalisatie A voltooid ===")
    println("Biomassa (vast)    : $FIXED_BIOMASSA cellen/L")
    println("MOI                : $(round(best_x[2],    digits=4))")
    println("Infectietijdstip   : $(round(best_x[1],    digits=2)) uur")
    println("Gewogen score      : $(round(best_score_nm, digits=6))")

    return (best_moi       = best_x[2],
            best_t_inf     = best_x[1],
            best_biomass   = FIXED_BIOMASSA,
            best_score     = best_score_nm)
end
