# ============================================================
#  CHAPTER 5: Metabolic cost and fixed-delay toxicity of Benzonase
#
#  This merges the original set of simulations (3a, 3b, 3e, 3f, 3g)
#  with the new toxicity investigation: all Parameters constructors
#  now use tau_death [h] (fixed delay between production start and
#  cell death) instead of the old continuous k_tox rate -- see the
#  updated Model2.jl / Model3.jl. The extended MOI sweep (3g) now
#  also tracks the lysogenic overtake time, and a new set of
#  representative low/medium/high MOI time series is added at the
#  end to directly check the overtake-time hypothesis.
# ============================================================
include("Model2.jl")
include("Model3.jl")
using Plots, Statistics, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

gr()

const BENZ_REF_3 = 0.001
const PN_REF_3   = 2.0e10

model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 1.32;  b = 170.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 2.26, 0.0, 10.0]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 1.71e-12;  duration = 40.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW  = 1.0e-12

tau_death = 1.0   # [h] -- fixed delay between production start and cell death (replaces k_tox)

naiveModel2   = Model2.loadFBAmodel(model_path)
all_ex_ids    = [id for id in keys(naiveModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2     = Model2.buildFbaCache(naiveModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2   = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                      "R_BIOMASS_Ec_iJO1366_core_53p95M")
naiveModel3   = Model3.loadFBAmodel(model_path)
lysogenModel3 = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
naiveFba3     = Model3.buildFbaCache(naiveModel3,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3   = Model3.buildFbaCache(lysogenModel3, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                      benz_id="R_BENZ_prod")

function make_p2(N0, t_inf, moi, f_prod_val=0.0015)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", tau_death, 0.0, 0.0,
        mu_max_vector, e_max_vector, f_prod_val, 0.05)
end

function make_p3(N0, t_inf, moi, f_prod_val=0.0015)
    Model3.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", tau_death, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod_val)
end

default_style = (tickfontsize=10, guidefontsize=11, legendfontsize=8,
                  bottom_margin=8Plots.mm, left_margin=10Plots.mm, top_margin=3Plots.mm)

# ------------------------------------------------------------
#  Safe axis-limit helpers: a plain MOI sweep can contain NaN entries
#  (failed/infeasible runs) or exact zeros (no lysogens/no benzonase at
#  very low MOI), both of which break automatic log-scale ranging or
#  leave a lot of dead white space on a linear axis. These filter such
#  values out before computing the plotted range.
# ------------------------------------------------------------
function safe_lin_ylims(vals; pad=1.15)
    finite_vals = filter(isfinite, vals)
    isempty(finite_vals) && return (0.0, 1.0)
    return (0.0, maximum(finite_vals) * pad)
end

function safe_log_ylims(vals; floor=1.0, pad=3.0)
    finite_vals = filter(x -> isfinite(x) && x > 0, vals)
    isempty(finite_vals) && return (floor, floor * 10)
    return (max(floor, minimum(finite_vals) / pad), maximum(finite_vals) * pad)
end

# ============================================================
#  3a. Lysogen growth rate vs f_prod (Model 3)
# ============================================================
println("=== 3a: Lysogen growth rate vs f_prod ===")

f_prod_values   = [0.0005, 0.001, 0.0015, 0.003, 0.005, 0.01]
mu_l_per_fp     = Float64[]
benz_eindwaarde = Float64[]

for fp in f_prod_values
    p_mu = make_p3(1e9, 999.0, 0.0, fp)
    sol_mu = Model3.run(p_mu)
    idx = argmin(abs.(sol_mu.t .- 3.0))
    u_stabiel = sol_mu.u[idx]
    Model3.fbaUpdate!(u_stabiel, p_mu)
    push!(mu_l_per_fp, p_mu.mu_l)

    sol_fp = Model3.run(make_p3(1e9, 2.0, 2.0, fp))
    if SciMLBase.successful_retcode(sol_fp)
        B_ts = [sol_fp.u[i][Model3.Benzind] for i in eachindex(sol_fp.u)]
        push!(benz_eindwaarde, maximum(B_ts))
    else
        push!(benz_eindwaarde, 0.0)
    end
end

pa = plot(f_prod_values, benz_eindwaarde, marker=:circle, lw=2, color=:green,
    xlabel="f_prod [-]", ylabel="Max Benzonase [mmol/L]", legend=false; default_style...)
pb = plot(f_prod_values, mu_l_per_fp, marker=:square, lw=2, color=:orange,
    xlabel="f_prod [-]", ylabel="Growth rate mu_l [1/h]", legend=false; default_style...)
fig3a = plot(pa, pb, layout=(1,2), size=(1000,440), margin=5Plots.mm)
savefig(fig3a, "h3a_fprod_effect.png")
println("  Saved: h3a_fprod_effect.png")

# ============================================================
#  3b. Benzonase time series: Model 2 vs Model 3
# ============================================================
println("\n=== 3b: Benzonase time series Model 2 vs Model 3 ===")

sol_m2 = Model2.run(make_p2(1e9, 2.0, 2.0, 0.0015))
sol_m3 = Model3.run(make_p3(1e9, 2.0, 2.0, 0.0015))

t_m2   = sol_m2.t; B_m2 = [sol_m2.u[i][Model2.Benzind] for i in eachindex(sol_m2.u)]
t_m3   = sol_m3.t; B_m3 = [sol_m3.u[i][Model3.Benzind] for i in eachindex(sol_m3.u)]
Lys_m3 = [sol_m3.u[i][Model3.lind] for i in eachindex(sol_m3.u)]
N_m3   = [sol_m3.u[i][Model3.Nind] for i in eachindex(sol_m3.u)]

p3b_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    legend=:outertop, legend_column=-1; default_style...)
plot!(p3b_benz, t_m2, B_m2, label="Growth-loss", color=:steelblue, lw=2, linestyle=:dash)
plot!(p3b_benz, t_m3, B_m3, label="FBA-cost",    color=:green,     lw=2)

p3b_pop = plot(xlabel="t [h]", ylabel="Cells/L",
    yscale=:log10, ylims=(1,:auto),
    legend=:outertop, legend_column=-1; default_style...)
plot!(p3b_pop, t_m3, max.(N_m3,1.0),   label="Naive",    color=:blue,  lw=2)
plot!(p3b_pop, t_m3, max.(Lys_m3,1.0), label="Lysogenic", color=:green, lw=2)

fig3b = plot(p3b_benz, p3b_pop, layout=(1,2), size=(1000,440), margin=5Plots.mm)
savefig(fig3b, "h3b_benzonase_tijdsreeks.png")
println("  Saved: h3b_benzonase_tijdsreeks.png")

# ============================================================
#  3e. MOI effect on Benzonase: Model 2 and Model 3
# ============================================================
println("\n=== 3e: MOI effect on Benzonase ===")

moi_sweep_benz  = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
benz_m2_per_moi = Float64[]; benz_m3_per_moi = Float64[]

for moi in moi_sweep_benz
    sol_m2_moi = Model2.run(make_p2(1e9, 2.0, moi, 0.0015))
    push!(benz_m2_per_moi, SciMLBase.successful_retcode(sol_m2_moi) ?
        maximum([sol_m2_moi.u[i][Model2.Benzind] for i in eachindex(sol_m2_moi.u)]) : 0.0)

    sol_m3_moi = Model3.run(make_p3(1e9, 2.0, moi, 0.0015))
    push!(benz_m3_per_moi, SciMLBase.successful_retcode(sol_m3_moi) ?
        maximum([sol_m3_moi.u[i][Model3.Benzind] for i in eachindex(sol_m3_moi.u)]) : 0.0)
end

fig3e = plot(xlabel="MOI [-]", ylabel="Max Benzonase [mmol/L]",
    xscale=:log10, legend=:outertop, legend_column=-1; default_style...)
plot!(fig3e, moi_sweep_benz, benz_m2_per_moi, marker=:circle, lw=2, color=:steelblue, label="Growth-loss")
plot!(fig3e, moi_sweep_benz, benz_m3_per_moi, marker=:square, lw=2, color=:green,     label="FBA-cost")
savefig(fig3e, "h3e_moi_effect_benzonase.png")
println("  Saved: h3e_moi_effect_benzonase.png")

# ============================================================
#  3f. Sensitivity analysis alpha_syn (Model 3)
# ============================================================
println("\n=== 3f: Sensitivity analysis alpha_syn ===")

alpha_syn_waarden = [0.5, 1.0, 2.0, 4.0]
colors_3f         = [:darkblue, :steelblue, :green, :darkorange]
labels_3f         = ["alpha=0.5", "alpha=1.0", "alpha=2.0 (ref)", "alpha=4.0"]

function make_p3_asyn(N0, t_inf, moi, asyn, f_prod_val=0.0015)
    # NOTE: e_max depends on alpha_syn and must be recomputed here, but
    # mu_max is a fixed biological property of each pathway and must stay
    # at mu_max_vector regardless of alpha_syn -- the earlier version of
    # this function passed the recomputed e_max in both slots by mistake.
    e_max_asyn = (asyn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
    Model3.Parameters(duration, N0, asyn, beta_deg, K_s, V_max, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", tau_death, 0.1, 0.0,
        mu_max_vector, e_max_asyn, f_prod_val)
end

p3f_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    legend=:outertop, legend_column=-1; default_style...)
p3f_lys  = plot(xlabel="t [h]", ylabel="Lysogenic cells/L",
    yscale=:log10, ylims=(1,:auto),
    legend=:outertop, legend_column=-1; default_style...)
p3f_emal = plot(xlabel="t [h]", ylabel="e_mal [-]",
    legend=:outertop, legend_column=-1; default_style...)

benz_asyn = Float64[]; t_max_asyn = Float64[]

for (idx, asyn) in enumerate(alpha_syn_waarden)
    sol_as = Model3.run(make_p3_asyn(1e9, 2.0, 2.0, asyn))
    if SciMLBase.successful_retcode(sol_as)
        t_as    = sol_as.t
        B_as    = [sol_as.u[i][Model3.Benzind] for i in eachindex(sol_as.u)]
        Lys_as  = max.([sol_as.u[i][Model3.lind]   for i in eachindex(sol_as.u)], 1.0)
        emal_as = [sol_as.u[i][Model3.Eind[2]]      for i in eachindex(sol_as.u)]
        i_max   = argmax(B_as)
        push!(benz_asyn, B_as[i_max]); push!(t_max_asyn, t_as[i_max])
        plot!(p3f_benz, t_as, B_as,    label=labels_3f[idx], color=colors_3f[idx], lw=2)
        plot!(p3f_lys,  t_as, Lys_as,  label=labels_3f[idx], color=colors_3f[idx], lw=2)
        plot!(p3f_emal, t_as, emal_as, label=labels_3f[idx], color=colors_3f[idx], lw=2)
    end
end

p3f_bar_benz = bar(string.(alpha_syn_waarden), benz_asyn,
    xlabel="alpha_syn [1/h]", ylabel="Max Benzonase [mmol/L]",
    color=colors_3f, legend=false; default_style...)
p3f_bar_t = bar(string.(alpha_syn_waarden), t_max_asyn,
    xlabel="alpha_syn [1/h]", ylabel="Time of max Benzonase [h]",
    color=colors_3f, legend=false; default_style...)

fig3f_ts  = plot(p3f_benz, p3f_lys, p3f_emal, layout=(1,3), size=(1650,480), margin=4Plots.mm)
fig3f_sum = plot(p3f_bar_benz, p3f_bar_t,     layout=(1,2), size=(1000,420), margin=5Plots.mm)
savefig(fig3f_ts,  "h3f_alphasyn_tijdsreeksen.png")
savefig(fig3f_sum, "h3f_alphasyn_samenvatting.png")
println("  Saved: h3f_alphasyn_tijdsreeksen.png | h3f_alphasyn_samenvatting.png")

# ============================================================
#  3g. Extended MOI sweep (Model 2 and Model 3), now also
#  tracking the lysogenic overtake time to test the toxicity
#  hypothesis directly (see header comment).
# ============================================================
println("\n=== 3g: Extended MOI sweep ===")

moi_sweep_g = [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.008,
               0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.5,
               1.0, 2.0, 5.0, 10.0]

benz_m2_g = Float64[]; lys_m2_g = Float64[]; biomass_loss_m2_g = Float64[]; tover_m2_g = Float64[]
benz_m3_g = Float64[]; lys_m3_g = Float64[]; biomass_loss_m3_g = Float64[]; tover_m3_g = Float64[]

function overtake_time(sol, N_idx, Lyt_idx, Lys_idx)
    N   = [sol.u[i][N_idx]   for i in eachindex(sol.u)]
    Lyt = [sol.u[i][Lyt_idx] for i in eachindex(sol.u)]
    Lys = [sol.u[i][Lys_idx] for i in eachindex(sol.u)]
    idx = findfirst(i -> Lys[i] > N[i] + Lyt[i], eachindex(sol.u))
    return isnothing(idx) ? NaN : sol.t[idx]
end

for moi in moi_sweep_g
    sol_m2_g = Model2.run(make_p2(1e9, 2.0, moi, 0.0015))
    if SciMLBase.successful_retcode(sol_m2_g)
        push!(benz_m2_g, maximum([sol_m2_g.u[i][Model2.Benzind] for i in eachindex(sol_m2_g.u)]))
        push!(lys_m2_g,  maximum([sol_m2_g.u[i][Model2.lind]    for i in eachindex(sol_m2_g.u)]))
        X_max_m2 = maximum([sol_m2_g.u[i][Model2.Nind] + sol_m2_g.u[i][Model2.Lind] +
                             sol_m2_g.u[i][Model2.lind] + sol_m2_g.u[i][Model2.Dind]
                             for i in eachindex(sol_m2_g.u)])
        X_end_m2 = sol_m2_g.u[end][Model2.Nind] + sol_m2_g.u[end][Model2.Lind] +
                   sol_m2_g.u[end][Model2.lind]  + sol_m2_g.u[end][Model2.Dind]
        push!(biomass_loss_m2_g, max(0.0, (X_max_m2 - X_end_m2) / X_max_m2))
        push!(tover_m2_g, overtake_time(sol_m2_g, Model2.Nind, Model2.Lind, Model2.lind))
    else
        push!(benz_m2_g, NaN); push!(lys_m2_g, NaN); push!(biomass_loss_m2_g, NaN); push!(tover_m2_g, NaN)
    end

    sol_m3_g = Model3.run(make_p3(1e9, 2.0, moi, 0.0015))
    if SciMLBase.successful_retcode(sol_m3_g)
        push!(benz_m3_g, maximum([sol_m3_g.u[i][Model3.Benzind] for i in eachindex(sol_m3_g.u)]))
        push!(lys_m3_g,  maximum([sol_m3_g.u[i][Model3.lind]    for i in eachindex(sol_m3_g.u)]))
        X_max_m3 = maximum([sol_m3_g.u[i][Model3.Nind] + sol_m3_g.u[i][Model3.Lind] +
                             sol_m3_g.u[i][Model3.lind] + sol_m3_g.u[i][Model3.Dind]
                             for i in eachindex(sol_m3_g.u)])
        X_end_m3 = sol_m3_g.u[end][Model3.Nind] + sol_m3_g.u[end][Model3.Lind] +
                   sol_m3_g.u[end][Model3.lind]  + sol_m3_g.u[end][Model3.Dind]
        push!(biomass_loss_m3_g, max(0.0, (X_max_m3 - X_end_m3) / X_max_m3))
        push!(tover_m3_g, overtake_time(sol_m3_g, Model3.Nind, Model3.Lind, Model3.lind))
    else
        push!(benz_m3_g, NaN); push!(lys_m3_g, NaN); push!(biomass_loss_m3_g, NaN); push!(tover_m3_g, NaN)
    end

    println("  MOI=$moi | M2: benz=$(round(benz_m2_g[end],sigdigits=3)) t_over=$(tover_m2_g[end]) | M3: benz=$(round(benz_m3_g[end],sigdigits=3)) t_over=$(tover_m3_g[end])")
end

function plot_lys_biomass(moi_vals, lys_vals, loss_vals, kleur)
    lys_ylims = safe_log_ylims(lys_vals)
    p = plot(moi_vals, lys_vals, xscale=:log10, yscale=:log10, ylims=lys_ylims,
        marker=:square, markersize=4, lw=2, color=kleur, label="Max lysogenic",
        xlabel="MOI [-]", ylabel="Max lysogenic cells/L",
        legend=:outertop, legend_column=-1; default_style...)
    plot!(twinx(p), moi_vals, loss_vals, xscale=:log10, ylims=(0.0, 1.05),
        marker=:diamond, markersize=4, lw=2, color=:tomato, linestyle=:dash,
        label="Biomass loss", ylabel="Relative biomass loss [-]", legend=false)
    return p
end

fig3g_m2 = plot(
    plot(moi_sweep_g, benz_m2_g, xscale=:log10, ylims=safe_lin_ylims(benz_m2_g),
        marker=:circle, lw=2, color=:steelblue,
        xlabel="MOI [-]", ylabel="Max Benzonase [mmol/L]", legend=false; default_style...),
    plot_lys_biomass(moi_sweep_g, lys_m2_g, biomass_loss_m2_g, :steelblue),
    layout=(1,2), size=(1100,460), margin=5Plots.mm)
savefig(fig3g_m2, "h3g_moi_sweep_M2.png")
println("  Saved: h3g_moi_sweep_M2.png")

fig3g_m3 = plot(
    plot(moi_sweep_g, benz_m3_g, xscale=:log10, ylims=safe_lin_ylims(benz_m3_g),
        marker=:circle, lw=2, color=:green,
        xlabel="MOI [-]", ylabel="Max Benzonase [mmol/L]", legend=false; default_style...),
    plot_lys_biomass(moi_sweep_g, lys_m3_g, biomass_loss_m3_g, :green),
    layout=(1,2), size=(1100,460), margin=5Plots.mm)
savefig(fig3g_m3, "h3g_moi_sweep_M3.png")
println("  Saved: h3g_moi_sweep_M3.png")

# t_overtake is still tracked and printed above for every MOI (see the
# println inside the sweep loop), it is just no longer plotted as a third
# panel here -- these two figures are back to the original 2-panel layout.

# ============================================================
#  NEW: representative time series at low / medium / high MOI,
#  to directly illustrate the overtake-time hypothesis:
#   - Low MOI: late overtake, substrate possibly already low -> little benz
#   - High MOI: early overtake, but fast collapse before peak biomass -> little benz
#   - An intermediate MOI is expected to allow both a reasonably early
#     takeover AND a reasonably large biomass at that point
# ============================================================
println("\n=== Representative time series (low / medium / high MOI) ===")
moi_examples = [0.001, 1.0, 10.0]
labels_ex    = ["Low MOI", "Medium MOI", "High MOI"]
colors_ex    = [:steelblue, :darkorange, :firebrick]

p_ex_bio  = plot(xlabel="t [h]", ylabel="Biomass [cells/L]",
    yscale=:log10, ylims=(1,:auto),
    legend=:outertop, legend_column=-1; default_style...)
p_ex_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    legend=:outertop, legend_column=-1; default_style...)

for (idx, moi) in enumerate(moi_examples)
    sol_ex = Model3.run(make_p3(1e9, 2.0, moi, 0.0015))
    if SciMLBase.successful_retcode(sol_ex)
        t_ex = sol_ex.t
        N_ex   = [sol_ex.u[i][Model3.Nind]   for i in eachindex(sol_ex.u)]
        Lyt_ex = [sol_ex.u[i][Model3.Lind]   for i in eachindex(sol_ex.u)]
        Lys_ex = [sol_ex.u[i][Model3.lind]   for i in eachindex(sol_ex.u)]
        B_ex   = [sol_ex.u[i][Model3.Benzind] for i in eachindex(sol_ex.u)]
        X_ex   = max.(N_ex .+ Lyt_ex .+ Lys_ex, 1.0)
        plot!(p_ex_bio,  t_ex, X_ex, label=labels_ex[idx], color=colors_ex[idx], lw=2)
        plot!(p_ex_benz, t_ex, B_ex, label=labels_ex[idx], color=colors_ex[idx], lw=2)
    end
end

fig_examples = plot(p_ex_bio, p_ex_benz, layout=(1,2), size=(1200,460), margin=5Plots.mm)
savefig(fig_examples, "h5_moi_examples_toxicity.png")
println("  Figure saved: h5_moi_examples_toxicity.png")

println("""

Compare the printed t_overtake / benz values from 3g (and the h5 example
figure) against the hypothesis: does the lowest MOI show a late t_overtake
combined with low benz? Does the highest MOI show an early t_overtake but
also a large biomass_loss and low benz? Is there a clear intermediate MOI
where benz is maximal? Report back what the numbers actually show so the
discussion section can be written around the real result rather than the
hypothesis.
""")

println("=== Chapter 5 simulations complete ===")