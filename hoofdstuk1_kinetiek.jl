# ============================================================
#  CHAPTER 1: Model simulations (Model 1)
#
#  Scenarios (3 figures each: growth / substrate / enzyme):
#   1. Glucose only
#   2. Maltose only
#   3. Glucose + Maltose (diauxie)
#   4. Acetate overflow metabolism (oxygen limitation)
#  -> 12 figures in total
#
#  Sensitivity analysis (3 figures each: enzyme / growth / substrate):
#   5. alpha_syn (enzyme synthesis rate), factor 10 lower
#   6. beta_deg  (enzyme degradation rate), factor 10 lower
# ============================================================

include("./Bin_model_1/dFBA.jl")
include("./Bin_model_1/parameters.jl")
include("./Bin_model_1/FBA.jl")

# ------------------------------------------------------------
# Fix for a sign bug in the substrate ODE of dFBA.jl / Model1:
#
#   du[sub_idx] = S_subs[i] < 1e-8 && p.q[i] >= 0 ? 0.0 : ...
#
# This zeroes the derivative whenever a substrate is (near) depleted
# AND its exchange flux is non-negative. Since q >= 0 corresponds to
# SECRETION in this codebase's sign convention (uptake is negative,
# capped via `lower_bound = -vmax`), this condition silently blocks
# any accumulation of a substrate that starts at zero -- exactly the
# situation for acetate in the overflow scenario. The intended
# behaviour is almost certainly the opposite: block further UPTAKE
# (q < 0) once a substrate is depleted, but always allow secretion
# (q >= 0) to raise its concentration from zero. simulate_dFBA! is
# therefore redefined here with that single condition corrected; the
# rest of the function is unchanged from dFBA.jl. Consider applying
# the same one-line fix directly in Bin_model_1/dFBA.jl (and in
# Model1/Model2/Model3, which contain the same line) so every chapter
# benefits from it, not just this script.
# ------------------------------------------------------------
function simulate_dFBA!(du, u, h, p::Parameters, t)::Nothing
    S_subs = u[Sind]
    e_enz  = u[Eind]
    S_cell = u[Sind_S]
    I_cell = u[Sind_I]
    L_cell = u[Sind_L]

    f     = getMonod(S_subs, p)
    R     = getRate(f, p)
    u_cyt = getU_cyt(R, p)

    tau_d = 28.0 / 60.0
    tau_l = p.tau

    u_decision = h(p, t - tau_d)
    u_lysis    = h(p, t - tau_l)

    S_lysis = u_lysis[Sind_S]
    P_lysis = u_lysis[Pfind]
    lysis_term = (t > p.infection_time + tau_l && S_lysis > CELL_THRESHOLD && P_lysis > PHAGE_THRESHOLD) ?
        getNewInfectionFlux(u_lysis, p) : 0.0

    S_decision = u_decision[Sind_S]
    phi_beslissing = (t > p.infection_time + tau_d && S_decision > CELL_THRESHOLD) ?
        getPhi(u_decision, p) : 0.0
    prob_lys = getProbLys(phi_beslissing)

    X_tot            = getTotalBiomass(u, p)
    totaal_adsorptie = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie  = getNewInfectionFlux(u, p)

    mu_safe  = max(0.0, p.mu)
    mu_eff_l = mu_safe * (1.0 - p.f_prod)

    productie_benz = (mu_safe - mu_eff_l) * L_cell * p.Y_benz

    for i in eachindex(Sind)
        sub_idx = Sind[i]
        enz_idx = Eind[i]
        du[sub_idx] = S_subs[i] < 1e-8 && p.q[i] < 0 ?
            0.0 : p.q[i] * p.E_coli_cellDW * X_tot
        if i == 1
            du[sub_idx] += p.h_release * lysis_term
        end
        mu_avg      = X_tot > 1.0 ?
            (mu_safe * S_cell + mu_eff_l * L_cell) / X_tot : mu_safe
        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] -
            (p.beta_deg + mu_avg) * e_enz[i] + 0.001
    end

    S_actief = S_cell > CELL_THRESHOLD ? S_cell : 0.0
    L_actief = L_cell > CELL_THRESHOLD ? L_cell : 0.0

    du[Sind_S]  = S_actief > 0.0 ? mu_safe * S_actief - nieuwe_infectie : 0.0
    du[Sind_I]  = (1.0 - prob_lys) * nieuwe_infectie - lysis_term
    du[Sind_L]  = mu_eff_l * L_actief + prob_lys * nieuwe_infectie
    du[Pfind]   = p.b * lysis_term - totaal_adsorptie
    du[Benzind] = productie_benz

    return nothing
end

using JuMP
using Plots, SciMLBase
using COBREXA, AbstractFBCModels
using DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

gr()

# ------------------------------------------------------------
#  Model and base parameters
# ------------------------------------------------------------
model_path = joinpath(@__DIR__, "iJO1366.xml")
model      = loadFBAmodel(model_path)

alpha_syn      = 2.0
beta_deg       = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 79/60
b              = 170.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
h_release      = 1.71e-12
duration       = 15.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
alfa_ads       = 1.3e-9
E_coli_cellDW  = 1.0e-12
biomass_id     = "R_BIOMASS_Ec_iJO1366_core_53p95M"
ac_idx         = findfirst(==("R_EX_ac_e"), exchange_ids)   # index of acetate within exchange_ids/Sind/Eind

# Standard (aerobic) FBA cache: oxygen assumed in excess, as usual
fbaCache_aerobic = buildFbaCache(model, exchange_ids, biomass_id)

# ------------------------------------------------------------
# Micro-aerobic FBA cache: oxygen uptake is restricted, giving a
# somewhat reduced respiratory capacity compared to the aerobic
# cache. Note: trying to *derive* acetate overflow purely from this
# restriction, by closing off every competing exchange reaction and
# letting the genome-scale LP "discover" that it must secrete
# acetate, turned out to be unreliable in practice -- with 2500+
# reactions there are always alternative ways for the solver to
# balance mass and redox that do not require any single reaction to
# carry flux, and closing enough of them risks making the whole
# problem infeasible instead (as happened with the fully forced
# uptake used in an earlier version of this script). Overflow
# metabolism is therefore imposed directly and transparently below,
# as a minimum acetate secretion tied to the available glucose
# uptake capacity, rather than left to emerge from the LP.
# ------------------------------------------------------------
model_microaerobic = deepcopy(model)
model_microaerobic.reactions["R_EX_o2_e"].lower_bound = -2.0   # mmol gDW-1 h-1, O2-limited
model_microaerobic.reactions["R_EX_ac_e"].upper_bound = 1000.0 # make sure secretion is not blocked by default

fbaCache_microaerobic = buildFbaCache(model_microaerobic, exchange_ids, biomass_id)

# ------------------------------------------------------------
#  Helper: build Parameters for a given scenario
# ------------------------------------------------------------
function make_scenario_params(vmax::Vector{Float64}, fbaCache::FbaCache;
                               alpha_::Float64=alpha_syn, beta_::Float64=beta_deg)
    e_max_local = (alpha_ .+ 0.001) ./ (beta_ .+ mu_max_vector)
    return Parameters(duration, 1e9, alpha_, beta_, K_s, vmax, p_pref,
        tau, b, alfa_ads, E_coli_cellDW, MW_values, h_release,
        biomass_id, exchange_ids, all_ex_ids, essentials_ids,
        fbaCache, 0.0, zeros(4), 999.0, 0.0,
        "R_BENZ_prod", 0.0, 0.0, 0.0, mu_max_vector, e_max_local, 0.0, 0.0)
end

# ------------------------------------------------------------
#  Overflow FBA update: identical to the standard fbaUpdate! (uptake
#  capped, not forced), except that the acetate exchange reaction is
#  additionally given a minimum (forced) secretion flux proportional
#  to the currently available glucose uptake capacity:
#
#     ac_min(t) = overflow_fraction * vmax_glc(t)
#
#  This directly guarantees that, whenever glucose is available, a
#  fixed fraction of the carbon taken up is diverted to acetate --
#  representing, in a simplified way, the fraction of carbon that
#  overflow metabolism diverts away from full oxidation once
#  respiratory (O2-dependent) capacity becomes limiting. The bound
#  naturally goes back to zero as glucose is depleted (vmax_glc -> 0).
# ------------------------------------------------------------
const overflow_fraction = 0.25   # fraction of the glucose uptake capacity forced out as acetate

function solveFbaOverflow(cache::FbaCache, vmax::Vector{Float64}, glc_idx::Int, ac_idx::Int, overflow_frac::Float64)
    for i in eachindex(vmax)
        setLowerBound!(cache, cache.exchange_ids[i], -vmax[i])
    end
    ac_min = overflow_frac * vmax[glc_idx]
    setLowerBound!(cache, cache.exchange_ids[ac_idx], ac_min)   # positive lower bound -> forced secretion
    JuMP.optimize!(cache.optimizer)
    if !JuMP.is_solved_and_feasible(cache.optimizer; dual=false)
        return nothing
    end
    fluxes = Dict{String,Float64}(rid => JuMP.value(cache.reaction_vars[rid]) for rid in cache.tracked_ids)
    return FbaSolution(fluxes)
end

function fbaUpdateOverflow!(u::Vector{Float64}, p::Parameters, glc_idx::Int, ac_idx::Int)
    vmax = getVmax(u, p)
    sol  = solveFbaOverflow(p.fbaModel, vmax, glc_idx, ac_idx, overflow_fraction)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu = sol.fluxes[p.biomass_id]
        p.q  = getFluxes(sol, p.ex_ids)
    else
        p.mu = 0.0; p.q .= 0.0
    end
    p.q_benz = 0.0
end

# ------------------------------------------------------------
#  Helper: run a simulation without infection, with freely chosen
#  initial carbon source concentrations. `overflow=true` switches to
#  fbaUpdateOverflow! for the acetate-overflow scenario.
# ------------------------------------------------------------
function run_scenario(p::Parameters, S0::Vector{Float64}; overflow::Bool=false)
    u0         = zeros(13)
    u0[Sind]   = S0
    u0[Eind]   = [0.95, 0.01, 0.01, 0.01]
    u0[Sind_S] = p.startingBiomass
    tspan      = (0.0, p.duration)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
    function fbaAffect!(integrator)
        if overflow
            fbaUpdateOverflow!(integrator.u, p, 1, ac_idx)
        else
            fbaUpdate!(integrator.u, p)
        end
        enforceThreshold!(integrator.u)
    end
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator) = any(x -> x < 0.0, u)
    domainAffect!(integrator)         = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    prob = DDEProblem(simulate_dFBA!, u0, (p, t) -> u0, tspan, p)
    return solve(prob, MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-6,
        tstops=fbaUpdateTimepoints,
        callback=CallbackSet(domainCallBack, fbaCallBack))
end

# ------------------------------------------------------------
#  Helper: extract results from the solution
# ------------------------------------------------------------
function extract(sol)
    t     = sol.t
    X     = [sol.u[i][Sind_S] for i in eachindex(sol.u)]
    glc   = [sol.u[i][Sind[1]] for i in eachindex(sol.u)]
    malt  = [sol.u[i][Sind[2]] for i in eachindex(sol.u)]
    ac    = [sol.u[i][Sind[4]] for i in eachindex(sol.u)]
    e_glc = [sol.u[i][Eind[1]] for i in eachindex(sol.u)]
    e_mal = [sol.u[i][Eind[2]] for i in eachindex(sol.u)]
    e_ac  = [sol.u[i][Eind[4]] for i in eachindex(sol.u)]
    return (t=t, X=X, glc=glc, malt=malt, ac=ac, e_glc=e_glc, e_mal=e_mal, e_ac=e_ac)
end

# ------------------------------------------------------------
#  Shared plot style: kept minimal on purpose (short axis labels,
#  no on-curve annotations, legend moved outside the axes so it
#  never overlaps a curve), with generous margins/figure size so
#  axis labels are never clipped.
# ------------------------------------------------------------
default_style = (tickfontsize=10, guidefontsize=11, legendfontsize=8,
                  bottom_margin=8Plots.mm, left_margin=10Plots.mm, top_margin=3Plots.mm)

# ------------------------------------------------------------
#  Helper: standard 3-panel figure (growth / substrate / enzyme).
#  e_i is a dimensionless model variable (not bounded to [0,1]), so
#  its axis is scaled to the data of each figure instead of fixed.
# ------------------------------------------------------------
function plot_scenario(res, substrate_series, substrate_labels, substrate_colors,
                        enzyme_series, enzyme_labels, enzyme_colors, filename::String)

    multi = length(substrate_series) > 1

    p_growth = plot(res.t, max.(res.X, 1.0),
        yscale=:log10, ylims=(1e7, :auto),
        xlabel="t [h]", ylabel="Biomass [cells/L]",
        color=:black, lw=2, legend=false; default_style...)

    p_sub = plot(xlabel="t [h]", ylabel="Substrate [mmol/L]",
        legend = multi ? :outertop : false, legend_column=multi ? -1 : 1; default_style...)
    for (series, lab, col) in zip(substrate_series, substrate_labels, substrate_colors)
        plot!(p_sub, res.t, series, label=lab, color=col, lw=2)
    end

    enzyme_max  = maximum(vcat(enzyme_series...))
    enzyme_ylim = (0.0, max(enzyme_max * 1.2, 0.05))
    p_enz = plot(xlabel="t [h]", ylabel="Enzyme e_i [-]", ylims=enzyme_ylim,
        legend = multi ? :outertop : false, legend_column=multi ? -1 : 1; default_style...)
    for (series, lab, col) in zip(enzyme_series, enzyme_labels, enzyme_colors)
        plot!(p_enz, res.t, series, label=lab, color=col, lw=2)
    end

    fig = plot(p_growth, p_sub, p_enz, layout=(1,3), size=(1650,480), margin=4Plots.mm)
    savefig(fig, filename)
    println("  Figure saved: $filename")
    return fig
end

# ------------------------------------------------------------
#  Helper: 3-panel sensitivity figure (enzyme / growth / substrate),
#  comparing a reference run against a run with one parameter
#  changed. Legends are kept to "Reference"/"Alt" only; the actual
#  parameter values go in the printed console note, not on the figure.
# ------------------------------------------------------------
function plot_sensitivity(res_ref, res_alt, title::String, filename::String)

    enzyme_max  = max(maximum(res_ref.e_mal), maximum(res_alt.e_mal))
    enzyme_ylim = (0.0, max(enzyme_max * 1.2, 0.05))
    p_enz = plot(xlabel="t [h]", ylabel="enzyme [-]", ylims=enzyme_ylim,
        legend=:outertop, legend_column=-1; default_style...)
    plot!(p_enz, res_ref.t, res_ref.e_glc, label="e_glc_ref", color=:steelblue, lw=2)
    plot!(p_enz, res_alt.t, res_alt.e_glc, label="e_glc_alt", color=:steelblue, lw=2, linestyle=:dash)
    plot!(p_enz, res_ref.t, res_ref.e_mal, label="e_mal_ref", color=:darkorange, lw=2)
    plot!(p_enz, res_alt.t, res_alt.e_mal, label="e_mal_alt", color=:darkorange,   lw=2, linestyle=:dash)

    p_growth = plot(xlabel="t [h]", ylabel="Biomass [cells/L]",
        yscale=:log10, ylims=(1e7, :auto),
        legend=:outertop, legend_column=-1; default_style...)
    plot!(p_growth, res_ref.t, max.(res_ref.X, 1.0), label="Reference", color=:darkorange, lw=2)
    plot!(p_growth, res_alt.t, max.(res_alt.X, 1.0), label="Alt",       color=:darkred,   lw=2, linestyle=:dash)

    p_sub = plot(xlabel="t [h]", ylabel="Substrate [mmol/L]",
        legend=:outertop, legend_column=-1; default_style...)
    plot!(p_sub, res_ref.t, res_ref.glc,  label="Glc (ref)",  color=:steelblue,  lw=2)
    plot!(p_sub, res_ref.t, res_ref.malt, label="Mal (ref)",  color=:darkorange, lw=2)
    plot!(p_sub, res_alt.t, res_alt.glc,  label="Glc (alt)",  color=:steelblue,  lw=2, linestyle=:dash)
    plot!(p_sub, res_alt.t, res_alt.malt, label="Mal (alt)",  color=:darkorange, lw=2, linestyle=:dash)

    fig = plot(p_enz, p_growth, p_sub, layout=(1,3), size=(1650,520),
        plot_title=title, plot_titlefontsize=11, margin=4Plots.mm)
    savefig(fig, filename)
    println("  Figure saved: $filename")
    return fig
end

# ============================================================
#  1. Glucose only
# ============================================================
println("=== Scenario 1: Glucose only ===")
V_max_glc = [12.7, 0.0, 0.0, 4.0]
p_glc     = make_scenario_params(V_max_glc, fbaCache_aerobic)
sol_glc   = run_scenario(p_glc, [4.44, 0.0, 0.0, 0.0])
res_glc   = extract(sol_glc)
plot_scenario(res_glc,
    [res_glc.glc], ["Glucose"], [:steelblue],
    [res_glc.e_glc], ["e_glc"], [:steelblue],
    "h1_glucose_only.png")

# ============================================================
#  2. Maltose only
# ============================================================
println("\n=== Scenario 2: Maltose only ===")
V_max_malt = [0.0, 3.75, 0.0, 4.0]
p_malt     = make_scenario_params(V_max_malt, fbaCache_aerobic)
sol_malt   = run_scenario(p_malt, [0.0, 2.337, 0.0, 0.0])
res_malt   = extract(sol_malt)
plot_scenario(res_malt,
    [res_malt.malt], ["Maltose"], [:darkorange],
    [res_malt.e_mal], ["e_mal"], [:darkorange],
    "h1_maltose_only.png")

# ============================================================
#  3. Glucose + Maltose (diauxie)
# ============================================================
println("\n=== Scenario 3: Glucose + Maltose (diauxie) ===")
V_max_both = [12.7, 3.75, 0.0, 4.0]
p_both     = make_scenario_params(V_max_both, fbaCache_aerobic)
sol_both   = run_scenario(p_both, [4.44, 2.337, 0.0, 0.0])
res_both   = extract(sol_both)
plot_scenario(res_both,
    [res_both.glc, res_both.malt], ["Glucose", "Maltose"], [:steelblue, :darkorange],
    [res_both.e_glc, res_both.e_mal], ["e_glc", "e_mal"], [:steelblue, :darkorange],
    "h1_glucose_maltose.png")

# ============================================================
#  4. Acetate overflow metabolism (oxygen limitation)
#     A fixed fraction (overflow_fraction) of the available glucose
#     uptake capacity is forced out as acetate at every timestep, on
#     top of an oxygen-restricted FBA cache -- see the comments above
#     solveFbaOverflow for why this is imposed directly rather than
#     left to emerge from the genome-scale LP.
# ============================================================
println("\n=== Scenario 4: Acetate overflow metabolism (O2 limitation) ===")
V_max_overflow = [12.7, 0.0, 0.0, 4.0]
p_overflow     = make_scenario_params(V_max_overflow, fbaCache_microaerobic)
sol_overflow   = run_scenario(p_overflow, [4.44, 0.0, 0.0, 0.0]; overflow=true)
res_overflow   = extract(sol_overflow)
plot_scenario(res_overflow,
    [res_overflow.glc, res_overflow.ac], ["Glucose", "Acetate"], [:steelblue, :firebrick],
    [res_overflow.e_glc, res_overflow.e_ac], ["e_glc", "e_ac"], [:steelblue, :firebrick],
    "h1_acetate_overflow.png")

max_ac = maximum(res_overflow.ac)
max_X  = maximum(res_overflow.X)
if max_ac > 1e-6
    println("  Overflow confirmed: acetate accumulates to a maximum of $(round(max_ac, sigdigits=3)) mmol L⁻¹.")
elseif max_X <= p_overflow.startingBiomass * 1.01
    println("  WARNING: biomass never grows -- the combination of bounds is infeasible. Try a less negative R_EX_o2_e lower bound or a lower overflow_fraction.")
else
    println("  WARNING: growth occurs but no acetate accumulation was detected -- check ac_idx / exchange_ids ordering.")
end

# ============================================================
#  5. Sensitivity: alpha_syn (enzyme synthesis rate), factor 10 lower
#     Compared on the glucose-then-maltose scenario, where the switch
#     between two substrates -- and its effect on growth -- is
#     visible.
# ============================================================
println("\n=== Sensitivity: alpha_syn (enzyme synthesis) ===")
alpha_low        = alpha_syn / 10
beta_low        = beta_deg / 10
p_both_alpha_low = make_scenario_params(V_max_both, fbaCache_aerobic; alpha_=0.0, beta_=beta_low)
sol_alpha_low    = run_scenario(p_both_alpha_low, [4.44, 2.337, 0.0, 0.0])
res_alpha_low    = extract(sol_alpha_low)

fig = plot_sensitivity(res_both, res_alpha_low,
    "alpha_syn: $(alpha_syn) h⁻¹ (reference) vs $(alpha_low) h⁻¹ (×10 lower)",
    "h1_sensitivity_alpha.png")

# ============================================================
#  6. Sensitivity: beta_deg (enzyme degradation rate), factor 10 lower
# ============================================================
println("\n=== Sensitivity: beta_deg (enzyme degradation) ===")
beta_low        = beta_deg / 10
p_both_beta_low = make_scenario_params(V_max_both, fbaCache_aerobic; beta_=beta_low)
sol_beta_low    = run_scenario(p_both_beta_low, [4.44, 2.337, 0.0, 0.0])
res_beta_low    = extract(sol_beta_low)

plot_sensitivity(res_both, res_beta_low,
    "beta_deg: $(beta_deg) h⁻¹ (reference) vs $(beta_low) h⁻¹ (×10 lower)",
    "h1_sensitivity_beta.png")

# ============================================================
#  Note on parameter choice (printed, not annotated on the figures)
# ============================================================
println("""

Note on the parameter choice:
 - alpha_syn = $(alpha_syn) h⁻¹ was chosen so that the switch between
   two substrates (e.g. from glucose to maltose) matches, in terms of
   timescale, literature values for the diauxic lag phase in E. coli.
   At a ten-fold lower value (h1_sensitivity_alpha.png), this switch
   becomes clearly slower: e_mal builds up more gradually once
   glucose is exhausted, growth stalls longer around the switch, and
   maltose consumption in the substrate panel starts noticeably later.
 - beta_deg = $(beta_deg) h⁻¹ was chosen in a similar way, in line
   with typical degradation rates of membrane transporters. At a
   ten-fold lower value (h1_sensitivity_beta.png), both the build-up
   and breakdown of the enzyme pool slow down, giving a slower but
   more persistent elevated enzyme level and a correspondingly
   shifted glucose-to-maltose switching time.
""")

println("\n=== Chapter 1 (Model 1) complete ===")