using Statistics
using SciMLBase
using Plots
using Printf
using JuMP
using HiGHS

# ============================================================
#  GEDEELDE OPTIMIZER UTILITIES
#  Te includen in alle model-specifieke optimizers via:
#    include("optimizer_utils.jl")
#
#  Bevat:
#   - bereken_benz_ref_fba(): correcte bovengrens via FBA
#   - nelder_mead():          generieke Nelder-Mead optimizer
#   - pareto_front():         Pareto-front over (Benzonase, P_max)
# ============================================================

# ============================================================
#  ANALYTISCHE BOVENGRENS BENZONASE VIA FBA
#
#  Berekent de maximale theoretische Benzonase-productie [mmol/L]
#  via één enkele FBA-solve waarbij:
#    - Biomassareactie = 0 (geen groei, alles naar Benzonase)
#    - Substraatopname = V_max (volledige capaciteit)
#    - Objective = maximaliseer R_BENZ_prod flux
#
#  Parameters:
#    benz_model  : FBA-model met R_BENZ_prod al toegevoegd
#                  (bv. lysogenModel3 na addBenzonase!)
#    exchange_ids: vector van exchange reactie IDs
#    V_max       : maximale opnamesnelheden [mmol/gDW/h]
#    biomass_id  : ID van de biomassareactie
#    benz_id     : ID van de Benzonase productiereactie
#    rho_DW      : celdrooggewicht [gDW/cel]
#    N0          : beginbiomassa [cellen/L]
#    duration    : simulatieduur [h]
#
#  Retourneert: BENZ_REF [mmol/L]
# ============================================================
function bereken_benz_ref_fba(benz_model,
                               exchange_ids::Vector{String},
                               V_max::Vector{Float64},
                               biomass_id::String,
                               benz_id::String,
                               rho_DW::Float64,
                               N0::Float64,
                               duration::Float64)

    # Bouw een tijdelijke FBA-cache met Benzonase als objective
    # Biomassareactie wordt op 0 gezet (geen groei)
    optimizer_ref = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(optimizer_ref)

    reaction_vars = Dict{String, JuMP.VariableRef}()
    for (rid, rxn) in benz_model.reactions
        lb = rid == biomass_id ? 0.0 : rxn.lower_bound
        ub = rid == biomass_id ? 0.0 : rxn.upper_bound
        reaction_vars[rid] = JuMP.@variable(optimizer_ref,
            lower_bound = lb,
            upper_bound = ub,
            base_name   = rid)
    end

    # Stoichiometrische balansen
    balances = Dict(mid => JuMP.AffExpr() for mid in keys(benz_model.metabolites))
    for (rid, rxn) in benz_model.reactions
        for (mid, coef) in rxn.stoichiometry
            JuMP.add_to_expression!(balances[mid], coef, reaction_vars[rid])
        end
    end
    for bal in values(balances)
        JuMP.@constraint(optimizer_ref, bal == 0.0)
    end

    # Substraatopname op V_max zetten (volledige capaciteit)
    for (i, eid) in enumerate(exchange_ids)
        if haskey(reaction_vars, eid)
            JuMP.set_lower_bound(reaction_vars[eid], -V_max[i])
        end
    end

    # Objective: maximaliseer Benzonase-productie
    JuMP.@objective(optimizer_ref, Max, reaction_vars[benz_id])
    JuMP.optimize!(optimizer_ref)

    if !JuMP.is_solved_and_feasible(optimizer_ref; dual=false)
        @warn "BENZ_REF FBA niet oplosbaar — gebruik veilige fallback van 1.0 mmol/L"
        return 1.0
    end

    # Flux in mmol/gDW/h → omzetten naar mmol/L over de volledige duratie
    q_benz_max = JuMP.value(reaction_vars[benz_id])   # mmol/gDW/h
    benz_ref   = q_benz_max * rho_DW * N0 * duration  # mmol/L

    println("  BENZ_REF (FBA, geen groei): $(round(benz_ref, sigdigits=4)) mmol/L")
    println("  q_benz_max = $(round(q_benz_max, sigdigits=4)) mmol/gDW/h | ",
            "N0=$(N0) cel/L | dur=$(duration) h")
    return benz_ref
end

# ============================================================
#  NELDER-MEAD SIMPLEX OPTIMIZER
#  Maximaliseert f (retourneert hogere score = beter)
#  Grenzen worden via clamping afgedwongen
# ============================================================
function nelder_mead(f, x0, delta, lb, ub;
                     max_iter=300, tol=1e-5,
                     alpha=1.0, gamma=2.0, rho=0.5, sigma=0.5)

    clamp_x(x) = [clamp(x[i], lb[i], ub[i]) for i in eachindex(x)]
    n = length(x0)

    simplex = Vector{Vector{Float64}}(undef, n+1)
    simplex[1] = clamp_x(x0)
    for i in 1:n
        v      = copy(x0)
        v[i]  += delta[i]
        simplex[i+1] = clamp_x(v)
    end
    scores = [f(v) for v in simplex]

    for iter in 1:max_iter
        ord     = sortperm(scores, rev=true)
        simplex = simplex[ord]
        scores  = scores[ord]

        if abs(scores[1] - scores[end]) < tol
            println("  Nelder-Mead geconvergeerd na $iter iteraties.")
            break
        end

        cog     = clamp_x(mean(simplex[1:n]))
        xr      = clamp_x(cog .+ alpha .* (cog .- simplex[end]))
        score_r = f(xr)

        if score_r > scores[1]
            xe      = clamp_x(cog .+ gamma .* (xr .- cog))
            score_e = f(xe)
            simplex[end] = score_e > score_r ? xe : xr
            scores[end]  = score_e > score_r ? score_e : score_r
        elseif score_r > scores[end]
            simplex[end] = xr
            scores[end]  = score_r
        else
            xc      = clamp_x(cog .+ rho .* (simplex[end] .- cog))
            score_c = f(xc)
            if score_c > scores[end]
                simplex[end] = xc
                scores[end]  = score_c
            else
                for i in 2:n+1
                    simplex[i] = clamp_x(simplex[1] .+ sigma .* (simplex[i] .- simplex[1]))
                    scores[i]  = f(simplex[i])
                end
            end
        end

        iter % 50 == 0 && println("  Iter $iter | beste score = $(round(scores[1], digits=6))")
    end

    return simplex[1], scores[1]
end

# ============================================================
#  PARETO-FRONT GENERATOR
#
#  Varieert het gewicht w_benz van 0.05 → 0.95.
#  Voor elk gewicht wordt één optimalisatierun gedraaid.
#  Daarna worden gedomineerde punten gefilterd.
#
#  Een punt A wordt gedomineerd door B als:
#    B.benz ≥ A.benz  EN  B.P_max ≤ A.P_max
#    (met minstens één strikte ongelijkheid)
#
#  Parameters:
#    run_single_opt : functie(w_benz::Float64) die een NamedTuple
#                    retourneert met velden:
#                      w, benz, P_max, t_inf, moi, N0, score
#    n_weights      : aantal gewichtspunten (default 11)
#    figname        : bestandsnaam voor de Pareto-plot
# ============================================================
function pareto_front(run_single_opt::Function;
                      n_weights::Int = 11,
                      figname::String = "pareto_front.png")

    weights = range(0.05, stop=0.95, length=n_weights)
    results = []

    println("\n=== Pareto-front berekening ($n_weights gewichten) ===")
    for (i, w) in enumerate(weights)
        println("\n--- Pareto punt $i/$n_weights | w_benz=$(round(w, digits=2)) ---")
        res = run_single_opt(w)
        push!(results, res)
        @printf("  → Benz=%.4f mmol/L | P_max=%.3e fagen/L | t_inf=%.2f h | MOI=%.4f | N0=%.2e\n",
            res.benz, res.P_max, res.t_inf, res.moi, res.N0)
    end

    # Filter niet-gedomineerde (Pareto-optimale) punten
    pareto_pts = filter(results) do pt_a
        !any(results) do pt_b
            pt_b.benz  >= pt_a.benz  &&
            pt_b.P_max <= pt_a.P_max &&
            (pt_b.benz > pt_a.benz || pt_b.P_max < pt_a.P_max)
        end
    end
    sort!(pareto_pts, by = x -> x.P_max)

    println("\n=== Pareto-front: $(length(pareto_pts)) niet-gedomineerde punten ===")
    println("  w_benz | Benzonase [mmol/L] | P_max [fagen/L] | t_inf [h] | MOI    | N0")
    println("  -------|-------------------|-----------------|-----------|--------|----------")
    for pt in pareto_pts
        @printf("  %-6.2f | %-17.5f | %-15.3e | %-9.2f | %-6.4f | %.2e\n",
            pt.w, pt.benz, pt.P_max, pt.t_inf, pt.moi, pt.N0)
    end

    # Plot alle gesimuleerde punten + Pareto-front
    all_benz = [r.benz  for r in results]
    all_pmax = [r.P_max for r in results]
    par_benz = [p.benz  for p in pareto_pts]
    par_pmax = [p.P_max for p in pareto_pts]

    fig = plot(all_benz, all_pmax,
        seriestype=:scatter,
        marker=:circle, markersize=5,
        color=:lightgray, label="Alle simulaties",
        xlabel="Max Benzonase [mmol L⁻¹]",
        ylabel="Max vrije fagen [fagen L⁻¹]",
        yscale=:log10,
        legend=:topright,
        bottom_margin=8Plots.mm, left_margin=10Plots.mm)

    plot!(fig, par_benz, par_pmax,
        seriestype=:scatter,
        marker=:star5, markersize=9,
        color=:steelblue, label="Pareto-front")
    plot!(fig, par_benz, par_pmax,
        lw=2, color=:steelblue, linestyle=:dash, label="")

    # Annoteer Pareto-punten met gewicht
    for pt in pareto_pts
        annotate!(fig, pt.benz, pt.P_max * 1.5,
            text("w=$(round(pt.w, digits=2))", 7, :steelblue, :center))
    end

    savefig(fig, figname)
    println("  Pareto-figuur opgeslagen: $figname")

    return pareto_pts
end
