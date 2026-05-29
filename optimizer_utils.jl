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
#   - pareto_grid_scan():     originele versie (ongewijzigd)
#   - pareto_grid_scan_fixed(): gecorrigeerde versie
# ============================================================

# ============================================================
#  ANALYTISCHE BOVENGRENS BENZONASE VIA FBA
# ============================================================
function bereken_benz_ref_fba(benz_model,
                               exchange_ids::Vector{String},
                               V_max::Vector{Float64},
                               biomass_id::String,
                               benz_id::String,
                               rho_DW::Float64,
                               N0::Float64,
                               duration::Float64)

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

    balances = Dict(mid => JuMP.AffExpr() for mid in keys(benz_model.metabolites))
    for (rid, rxn) in benz_model.reactions
        for (mid, coef) in rxn.stoichiometry
            JuMP.add_to_expression!(balances[mid], coef, reaction_vars[rid])
        end
    end
    for bal in values(balances)
        JuMP.@constraint(optimizer_ref, bal == 0.0)
    end

    for (i, eid) in enumerate(exchange_ids)
        if haskey(reaction_vars, eid)
            JuMP.set_lower_bound(reaction_vars[eid], -V_max[i])
        end
    end

    JuMP.@objective(optimizer_ref, Max, reaction_vars[benz_id])
    JuMP.optimize!(optimizer_ref)

    if !JuMP.is_solved_and_feasible(optimizer_ref; dual=false)
        @warn "BENZ_REF FBA niet oplosbaar — gebruik veilige fallback van 1.0 mmol/L"
        return 1.0
    end

    q_benz_max = JuMP.value(reaction_vars[benz_id])
    benz_ref   = q_benz_max * rho_DW * N0 * duration

    println("  BENZ_REF (FBA, geen groei): $(round(benz_ref, sigdigits=4)) mmol/L")
    println("  q_benz_max = $(round(q_benz_max, sigdigits=4)) mmol/gDW/h | ",
            "N0=$(N0) cel/L | dur=$(duration) h")
    return benz_ref
end

# ============================================================
#  NELDER-MEAD SIMPLEX OPTIMIZER
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
#  ORIGINELE PARETO GRID SCAN (ongewijzigd)
# ============================================================
function pareto_grid_scan(p_initial, run_func::Function,
                           Benzind, Pfind,
                           moi_values, t_inf_values,
                           maak_parameters::Function;
                           figbase::String = "pareto_grid")

    results = []
    total   = length(moi_values) * length(t_inf_values)
    counter = 0

    println("\n=== Pareto grid scan ($total simulaties) ===")
    for moi in moi_values
        for t_inf in t_inf_values
            counter += 1
            p_test = maak_parameters(p_initial, moi, t_inf)
            try
                sol = redirect_stdout(devnull) do
                    redirect_stderr(devnull) do
                        run_func(p_test)
                    end
                end
                !SciMLBase.successful_retcode(sol) && continue
                max_benz     = maximum(sol[Benzind, :])
                P_max        = maximum(sol[Pfind, :])
                infectiedosis = moi * p_initial.startingBiomass
                idx_95       = findfirst(sol[Benzind, :] .>= 0.95 * max_benz)
                t_95         = isnothing(idx_95) ? p_initial.duration : sol.t[idx_95]
                @printf("  [%3d/%d] MOI=%-6.3f | t_inf=%-5.1f | benz=%.4e | P_max=%.3e | t95=%.2f\n",
                    counter, total, moi, t_inf, max_benz, P_max, t_95)
                push!(results, (moi=moi, t_inf=t_inf, max_benz=max_benz,
                                P_max=P_max, infectiedosis=infectiedosis, t_95=t_95))
            catch e
                @warn "Simulatie gefaald moi=$moi, t_inf=$t_inf: $e"
            end
        end
    end

    valid = filter(r -> isfinite(r.max_benz) && r.max_benz > 0, results)

    all_benz  = [r.max_benz      for r in valid]
    all_dose  = [r.infectiedosis  for r in valid]
    all_t95   = [r.t_95          for r in valid]
    all_tinf  = [r.t_inf         for r in valid]
    all_moi   = [r.moi           for r in valid]

    function get_pareto_front(xs, ys, maximize_x, minimize_y)
        pts = collect(zip(xs, ys, 1:length(xs)))
        pareto = filter(pts) do (ax, ay, _)
            !any(pts) do (bx, by, _)
                (maximize_x ? bx >= ax : bx <= ax) &&
                (minimize_y ? by <= ay : by >= ay) &&
                ((maximize_x ? bx > ax : bx < ax) ||
                 (minimize_y ? by < ay : by > ay))
            end
        end
        sort(pareto, by = p -> p[1])
    end

    par1 = get_pareto_front(all_benz, all_dose, true, true)
    px1  = [p[1] for p in par1]
    py1  = [p[2] for p in par1]

    fig1 = scatter(all_benz, all_dose,
        zcolor         = all_tinf,
        marker         = :circle, markersize = 6,
        colorbar_title = "t_inf [h]",
        label          = "All simulations",
        xlabel         = "Max Benzonase [mmol L⁻¹]",
        ylabel         = "Infection dose [phages L⁻¹]",
        legend         = :topright,
        bottom_margin  = 8Plots.mm, left_margin = 10Plots.mm)
    plot!(fig1, px1, py1, lw=2, color=:red, linestyle=:dash, label="Pareto front")
    scatter!(fig1, px1, py1, marker=:star5, markersize=9, color=:red, label="Pareto points")
    savefig(fig1, figbase * "_benz_vs_dosis.png")

    par2 = get_pareto_front(all_benz, all_t95, true, true)
    px2  = [p[1] for p in par2]
    py2  = [p[2] for p in par2]

    fig2 = scatter(all_benz, all_t95,
        zcolor         = all_moi,
        marker         = :circle, markersize = 6,
        colorbar_title = "MOI",
        label          = "All simulations",
        xlabel         = "Max Benzonase [mmol L⁻¹]",
        ylabel         = "t₉₅ [h]",
        legend         = :topright,
        bottom_margin  = 8Plots.mm, left_margin = 10Plots.mm)
    plot!(fig2, px2, py2, lw=2, color=:red, linestyle=:dash, label="Pareto front")
    scatter!(fig2, px2, py2, marker=:star5, markersize=9, color=:red, label="Pareto points")
    savefig(fig2, figbase * "_benz_vs_t95.png")

    return valid
end

# ============================================================
#  GECORRIGEERDE PARETO GRID SCAN
#
#  Wijzigingen t.o.v. origineel:
#   - Front 1: y-as = (t₉₅ - t_inf), de netto productietijd
#              gemeten vanaf infectie → echte procesvariabele
#   - Front 2: y-as = MOI (input als kostenmaatstaf),
#              want minimale infectiedosis was de doelstelling
#              in de optimizer, niet P_max
#   - Geen kleurencode: punten in eenvoudige grijstint,
#              de front-lijn en drie gelabelde punten
#              dragen alle relevante informatie
#   - Geen titels op de grafieken
#   - Alle labels in het Engels
#   - Drie representatieve punten aangeduid op elk front:
#              A = maximale Benzonase (rechts op front)
#              B = minimale y-waarde (links op front)
#              C = compromispunt (midden van front)
# ============================================================
function pareto_grid_scan_fixed(
        p_initial,
        run_func::Function,
        Benzind::Int,
        Pfind::Int,
        Nind::Int,
        moi_values::Vector{Float64},
        t_inf_values::Vector{Float64},
        maak_parameters::Function;
        figbase::String = "pareto_grid_fixed")

    results = []
    total   = length(moi_values) * length(t_inf_values)
    counter = 0

    println("\n=== Pareto grid scan ($total simulations) ===")
    println("    Front 1: Max Benzonase  vs  (t₉₅ - t_inf)  [production vs net process time]")
    println("    Front 2: Max Benzonase  vs  MOI             [production vs infection dose]")

    for moi in moi_values
        for t_inf in t_inf_values
            counter += 1
            p_test = maak_parameters(p_initial, moi, t_inf)
            try
                sol = redirect_stdout(devnull) do
                    redirect_stderr(devnull) do
                        run_func(p_test)
                    end
                end
                !SciMLBase.successful_retcode(sol) && continue

                B_ts     = sol[Benzind, :]
                max_benz = maximum(B_ts)
                max_benz < 1e-8 && continue

                # Netto productietijd: tijd vanaf infectie tot 95% van max Benzonase
                # Dit is een echte procesvariabele die geminimaliseerd wordt
                idx_95      = findfirst(B_ts .>= 0.95 * max_benz)
                t_95        = isnothing(idx_95) ? p_initial.duration : sol.t[idx_95]
                dt_prod     = max(t_95 - t_inf, 0.0)

                @printf("  [%3d/%d] MOI=%-6.3f | t_inf=%-4.1f | Benz=%.5f | dt_prod=%5.2f h\n",
                    counter, total, moi, t_inf, max_benz, dt_prod)

                push!(results, (
                    moi      = moi,
                    t_inf    = t_inf,
                    max_benz = max_benz,
                    dt_prod  = dt_prod,
                ))
            catch e
                @warn "Simulation failed moi=$moi, t_inf=$t_inf: $e"
            end
        end
    end

    isempty(results) && (@warn "No valid results!"; return results)

    # ----------------------------------------------------------
    #  Pareto-dominantie filter
    #  Punt A domineert punt B als:
    #    A.benz >= B.benz  (maximaliseren)
    #    A.obj2 <= B.obj2  (minimaliseren)
    #  met minstens één strikte ongelijkheid.
    # ----------------------------------------------------------
    function is_pareto(xs_max::Vector{Float64}, ys_min::Vector{Float64})
        n = length(xs_max)
        dominated = fill(false, n)
        for i in 1:n
            for j in 1:n
                i == j && continue
                if xs_max[j] >= xs_max[i] && ys_min[j] <= ys_min[i] &&
                   (xs_max[j] > xs_max[i]  || ys_min[j] < ys_min[i])
                    dominated[i] = true
                    break
                end
            end
        end
        return .!dominated
    end

    all_benz    = [r.max_benz for r in results]
    all_dtprod  = [r.dt_prod  for r in results]
    all_moi     = [r.moi      for r in results]
    all_tinf    = [r.t_inf    for r in results]

    # ----------------------------------------------------------
    #  Drie representatieve punten selecteren op een front
    #  A = maximale Benzonase (rechts)
    #  B = minimale y-waarde (links/onder)
    #  C = compromispunt (midden van de front-indices)
    # ----------------------------------------------------------
    function select_three(mask, xs, ys)
        idx   = findall(mask)
        ord   = sortperm(xs[idx])
        idx_s = idx[ord]
        n     = length(idx_s)
        A = idx_s[end]              # maximale Benzonase
        B = idx_s[argmin(ys[idx_s])]# minimale y
        C = idx_s[max(1, n ÷ 2)]   # midden
        return A, B, C
    end

    # ----------------------------------------------------------
    #  FIGUUR 1: Benzonase vs netto productietijd (t₉₅ - t_inf)
    #
    #  Trade-off: hogere Benzonase-productie vereist meer
    #  lysogenen, die langer nodig hebben om op te bouwen.
    #  Vroege infectie bij lage MOI geeft meer lysogenen maar
    #  een langere netto productietijd. Hoge MOI geeft snelle
    #  lysis maar een lager productieniveau.
    # ----------------------------------------------------------
    par1_mask = is_pareto(all_benz, all_dtprod)
    par1_ord  = sortperm(all_benz[par1_mask])
    par1_x    = all_benz[par1_mask][par1_ord]
    par1_y    = all_dtprod[par1_mask][par1_ord]

    A1, B1, C1 = select_three(par1_mask, all_benz, all_dtprod)

    fig1 = scatter(all_benz, all_dtprod,
        marker            = :circle,
        markersize        = 5,
        markercolor       = :gray60,
        markerstrokewidth = 0,
        label             = "All simulations",
        xlabel            = "Max Benzonase [mmol L⁻¹]",
        ylabel            = "t₉₅ - t_inf [h]",
        legend            = :topright,
        bottom_margin     = 8Plots.mm,
        left_margin       = 10Plots.mm,
        top_margin        = 6Plots.mm,
        size              = (620, 500))

    plot!(fig1, par1_x, par1_y,
        lw        = 2.5,
        color     = :red,
        linestyle = :dash,
        label     = "Pareto front")

    scatter!(fig1, all_benz[par1_mask], all_dtprod[par1_mask],
        marker            = :circle,
        markersize        = 7,
        markercolor       = :red,
        markerstrokewidth = 0,
        label             = "Pareto points")

    # Drie punten annoteren
    for (lbl, idx) in [("A", A1), ("B", B1), ("C", C1)]
        scatter!(fig1, [all_benz[idx]], [all_dtprod[idx]],
            marker      = :star5,
            markersize  = 12,
            markercolor = :darkred,
            markerstrokewidth = 0.5,
            label       = "")
        annotate!(fig1, all_benz[idx], all_dtprod[idx] + 0.15,
            text(lbl, 10, :darkred, :center))
    end

    savefig(fig1, figbase * "_benz_vs_dtprod.png")
    println("\nFigure saved: $(figbase)_benz_vs_dtprod.png")

    # ----------------------------------------------------------
    #  FIGUUR 2: Benzonase vs MOI
    #
    #  Trade-off: hogere MOI leidt tot meer directe infecties
    #  en potentieel meer Benzonase via lytische route, maar
    #  de infectiedosis is een directe procesinput die minimaal
    #  gehouden moet worden (kosten, faagreiniging).
    # ----------------------------------------------------------
    par2_mask = is_pareto(all_benz, all_moi)
    par2_ord  = sortperm(all_benz[par2_mask])
    par2_x    = all_benz[par2_mask][par2_ord]
    par2_y    = all_moi[par2_mask][par2_ord]

    A2, B2, C2 = select_three(par2_mask, all_benz, all_moi)

    fig2 = scatter(all_benz, all_moi,
        marker            = :circle,
        markersize        = 5,
        markercolor       = :gray60,
        markerstrokewidth = 0,
        label             = "All simulations",
        xlabel            = "Max Benzonase [mmol L⁻¹]",
        ylabel            = "MOI [-]",
        legend            = :topright,
        bottom_margin     = 8Plots.mm,
        left_margin       = 10Plots.mm,
        top_margin        = 6Plots.mm,
        size              = (620, 500))

    plot!(fig2, par2_x, par2_y,
        lw        = 2.5,
        color     = :red,
        linestyle = :dash,
        label     = "Pareto front")

    scatter!(fig2, all_benz[par2_mask], all_moi[par2_mask],
        marker            = :circle,
        markersize        = 7,
        markercolor       = :red,
        markerstrokewidth = 0,
        label             = "Pareto points")

    for (lbl, idx) in [("A", A2), ("B", B2), ("C", C2)]
        scatter!(fig2, [all_benz[idx]], [all_moi[idx]],
            marker      = :star5,
            markersize  = 12,
            markercolor = :darkred,
            markerstrokewidth = 0.5,
            label       = "")
        annotate!(fig2, all_benz[idx], all_moi[idx] * 1.08,
            text(lbl, 10, :darkred, :center))
    end

    savefig(fig2, figbase * "_benz_vs_moi.png")
    println("Figure saved: $(figbase)_benz_vs_moi.png")

    # ----------------------------------------------------------
    #  Gecombineerde figuur
    # ----------------------------------------------------------
    fig_combined = plot(fig1, fig2,
        layout = (1, 2),
        size   = (1250, 500),
        margin = 8Plots.mm)
    savefig(fig_combined, figbase * "_combined.png")
    println("Figure saved: $(figbase)_combined.png")

    # ----------------------------------------------------------
    #  Tekstueel overzicht van de drie punten per front
    # ----------------------------------------------------------
    println("\n=== Pareto points: Benzonase vs net process time ===")
    println("  Label | MOI       | t_inf [h] | Benzonase [mmol/L] | dt_prod [h]")
    for (lbl, idx) in [("A", A1), ("C", C1), ("B", B1)]
        @printf("  %-5s | %-9.4f | %-9.1f | %-18.5f | %.2f\n",
            lbl, results[idx].moi, results[idx].t_inf,
            results[idx].max_benz, results[idx].dt_prod)
    end

    println("\n=== Pareto points: Benzonase vs MOI ===")
    println("  Label | MOI       | t_inf [h] | Benzonase [mmol/L]")
    for (lbl, idx) in [("A", A2), ("C", C2), ("B", B2)]
        @printf("  %-5s | %-9.4f | %-9.1f | %.5f\n",
            lbl, results[idx].moi, results[idx].t_inf,
            results[idx].max_benz)
    end

    return results
end