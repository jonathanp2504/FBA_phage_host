# ============================================================
#  HOOFDSTUK 3: Metabolic cost of Benzonase
# ============================================================
include("Model2.jl")
include("Model3.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

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
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.35, 0.0, 0.0,
        mu_max_vector, e_max_vector, f_prod_val, 0.05)
end

function make_p3(N0, t_inf, moi, f_prod_val=0.0015)
    Model3.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.35, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod_val)
end

# ============================================================
#  3a. Lysogen growth rate vs f_prod (Model 3)
# ============================================================
println("=== 3a: Lysogen growth rate vs f_prod ===")

f_prod_values   = [0.0005, 0.001, 0.0015, 0.003, 0.005, 0.01]
mu_l_per_fp     = Float64[]
benz_eindwaarde = Float64[]

for fp in f_prod_values
    # Exacte mu_l via korte simulatie zonder infectie op stabiel groeipunt
    p_mu = make_p3(1e9, 999.0, 0.0, fp)
    sol_mu = Model3.run(p_mu)
    idx = argmin(abs.(sol_mu.t .- 3.0))
    u_stabiel = sol_mu.u[idx]
    Model3.fbaUpdate!(u_stabiel, p_mu)
    push!(mu_l_per_fp, p_mu.mu_l)

    # Benzonase via infectiesimulatie
    sol_fp = Model3.run(make_p3(1e9, 2.0, 2.0, fp))
    if SciMLBase.successful_retcode(sol_fp)
        B_ts = [sol_fp.u[i][Model3.Benzind] for i in eachindex(sol_fp.u)]
        push!(benz_eindwaarde, maximum(B_ts))
    else
        push!(benz_eindwaarde, 0.0)
    end
end

pa = plot(f_prod_values, benz_eindwaarde,
    label="Max Benzonase (Model 3)",
    marker=:circle, lw=2, color=:green,
    xlabel="f_prod [-]", ylabel="Max Benzonase [mmol L⁻¹]",
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
pb = plot(f_prod_values, mu_l_per_fp,
    label="Growth rate µ_l (Model 3)",
    marker=:square, lw=2, color=:orange,
    xlabel="f_prod [-]", ylabel="Growth rate µ_l [h⁻¹]",
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
fig3a = plot(pa, pb, layout=(1,2), size=(900,420), margin=6Plots.mm)
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

p3b_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p3b_benz, t_m2, B_m2, label="Model 2 (growth-loss)", color=:steelblue, lw=2, linestyle=:dash)
plot!(p3b_benz, t_m3, B_m3, label="Model 3 (FBA)",          color=:green,     lw=2)

p3b_pop = plot(xlabel="t [h]", ylabel="Cells L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p3b_pop, t_m3, [N_m3 Lys_m3], label=["Naive N" "Lysogen l"],
    color=[:blue :green], lw=2)

fig3b = plot(p3b_benz, p3b_pop, layout=(1,2), size=(900,430), margin=6Plots.mm)
savefig(fig3b, "h3b_benzonase_tijdsreeks.png")
println("  Saved: h3b_benzonase_tijdsreeks.png")

sol_m3 = Model3.run(make_p3(1e9, 2.0, 2.0, 0.0015))

# ============================================================
#  3e. MOI effect on Benzonase: Model 2 and Model 3
# ============================================================
println("\n=== 3e: MOI effect on Benzonase ===")

moi_sweep_benz = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
benz_m2_per_moi = Float64[]; benz_m3_per_moi = Float64[]

for moi in moi_sweep_benz
    sol_m2_moi = Model2.run(make_p2(1e9, 2.0, moi, 0.0015))
    push!(benz_m2_per_moi, SciMLBase.successful_retcode(sol_m2_moi) ?
        maximum([sol_m2_moi.u[i][Model2.Benzind] for i in eachindex(sol_m2_moi.u)]) : 0.0)

    sol_m3_moi = Model3.run(make_p3(1e9, 2.0, moi, 0.0015))
    push!(benz_m3_per_moi, SciMLBase.successful_retcode(sol_m3_moi) ?
        maximum([sol_m3_moi.u[i][Model3.Benzind] for i in eachindex(sol_m3_moi.u)]) : 0.0)
end

fig3e = plot(xlabel="MOI [-]", ylabel="Max Benzonase [mmol L⁻¹]",
    xscale=:log10, legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
plot!(fig3e, moi_sweep_benz, benz_m2_per_moi,
    marker=:circle, lw=2, color=:steelblue, label="Model 2 (growth-loss)")
plot!(fig3e, moi_sweep_benz, benz_m3_per_moi,
    marker=:square, lw=2, color=:green,     label="Model 3 (FBA)")
savefig(fig3e, "h3e_moi_effect_benzonase.png")
println("  Saved: h3e_moi_effect_benzonase.png")

# ============================================================
#  3f. Sensitivity analysis alpha_syn (Model 3)
# ============================================================
println("\n=== 3f: Sensitivity analysis alpha_syn ===")

alpha_syn_waarden = [0.5, 1.0, 2.0, 4.0]
kleuren_3f        = [:darkblue :steelblue :green :darkorange]
labels_3f         = ["α_syn=0.5" "α_syn=1.0" "α_syn=2.0 (ref)" "α_syn=4.0"]

function make_p3_asyn(N0, t_inf, moi, asyn, f_prod_val=0.0015)
    e_max_asyn = (asyn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
    Model3.Parameters(duration, N0, asyn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.05, 0.1, 0.0,
        e_max_asyn, e_max_asyn, f_prod_val)
end

p3f_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
p3f_lys  = plot(xlabel="t [h]", ylabel="Lysogenic cells L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
p3f_emal = plot(xlabel="t [h]", ylabel="e_mal (LamB proxy) [-]",
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)

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
        plot!(p3f_benz, t_as, B_as,    label=labels_3f[idx], color=kleuren_3f[idx], lw=2)
        plot!(p3f_lys,  t_as, Lys_as,  label=labels_3f[idx], color=kleuren_3f[idx], lw=2)
        plot!(p3f_emal, t_as, emal_as, label=labels_3f[idx], color=kleuren_3f[idx], lw=2)
    end
end

p3f_bar_benz = bar(string.(alpha_syn_waarden), benz_asyn,
    xlabel="α_syn [h⁻¹]", ylabel="Max Benzonase [mmol L⁻¹]",
    color=collect(kleuren_3f), legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
p3f_bar_t = bar(string.(alpha_syn_waarden), t_max_asyn,
    xlabel="α_syn [h⁻¹]", ylabel="Time of max Benzonase [h]",
    color=collect(kleuren_3f), legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)

fig3f_ts  = plot(p3f_benz, p3f_lys, p3f_emal, layout=(1,3), size=(1400,430), margin=6Plots.mm)
fig3f_sum = plot(p3f_bar_benz, p3f_bar_t,      layout=(1,2), size=(900,400),  margin=6Plots.mm)
savefig(fig3f_ts,  "h3f_alphasyn_tijdsreeksen.png")
savefig(fig3f_sum, "h3f_alphasyn_samenvatting.png")
println("  Saved: h3f_alphasyn_tijdsreeksen.png | h3f_alphasyn_samenvatting.png")

# ============================================================
#  3g UITGEBREID: meer MOI punten
# ============================================================
moi_sweep_g  = [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.008,
                0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.5,
                1.0, 2.0, 5.0, 10.0]
benz_m2_g    = Float64[]; lys_m2_g = Float64[]; biomass_loss_m2_g = Float64[]
benz_m3_g    = Float64[]; lys_m3_g = Float64[]; biomass_loss_m3_g = Float64[]

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
    else
        push!(benz_m2_g, NaN); push!(lys_m2_g, NaN); push!(biomass_loss_m2_g, NaN)
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
    else
        push!(benz_m3_g, NaN); push!(lys_m3_g, NaN); push!(biomass_loss_m3_g, NaN)
    end
end

function plot_lys_biomass(moi_vals, lys_vals, loss_vals, kleur)
    p = plot(moi_vals, lys_vals,
        xscale        = :log10,
        yscale        = :log10,
        marker        = :square,
        markersize    = 4,
        linewidth     = 2,
        color         = kleur,
        label         = "Max lysogenic cells",
        xlabel        = "Initial MOI [-]",
        ylabel        = "Max lysogenic cells [cells L⁻¹]",
        tickfontsize  = 11,
        guidefontsize = 13,
        legendfontsize= 10,
        legend        = :topleft,
        bottom_margin = 8Plots.mm,
        left_margin   = 10Plots.mm,
        right_margin  = 14Plots.mm)

    plot!(twinx(p), moi_vals, loss_vals,
        xscale        = :log10,
        marker        = :diamond,
        markersize    = 4,
        linewidth     = 2,
        color         = :tomato,
        linestyle     = :dash,
        label         = "Relative biomass loss",
        ylabel        = "Relative biomass loss [-]",
        tickfontsize  = 11,
        guidefontsize = 13,
        legendfontsize= 10,
        legend        = :bottomright)

    return p
end

fig3g_m2 = plot(
    plot(moi_sweep_g, benz_m2_g,
        xscale        = :log10,
        marker        = :circle,
        markersize    = 4,
        linewidth     = 2,
        color         = :steelblue,
        xlabel        = "Initial MOI [-]",
        ylabel        = "Max Benzonase [mmol L⁻¹]",
        tickfontsize  = 11,
        guidefontsize = 13,
        legend        = false,
        bottom_margin = 8Plots.mm,
        left_margin   = 10Plots.mm),
    plot_lys_biomass(moi_sweep_g, lys_m2_g, biomass_loss_m2_g, :steelblue),
    layout=(1,2), size=(1000,420))
savefig(fig3g_m2, "h3g_moi_sweep_M2.png")
println("  Saved: h3g_moi_sweep_M2.png")

fig3g_m3 = plot(
    plot(moi_sweep_g, benz_m3_g,
        xscale        = :log10,
        marker        = :circle,
        markersize    = 4,
        linewidth     = 2,
        color         = :green,
        xlabel        = "Initial MOI [-]",
        ylabel        = "Max Benzonase [mmol L⁻¹]",
        tickfontsize  = 11,
        guidefontsize = 13,
        legend        = false,
        bottom_margin = 8Plots.mm,
        left_margin   = 10Plots.mm),
    plot_lys_biomass(moi_sweep_g, lys_m3_g, biomass_loss_m3_g, :green),
    layout=(1,2), size=(1000,420))
savefig(fig3g_m3, "h3g_moi_sweep_M3.png")
println("  Saved: h3g_moi_sweep_M3.png")

println("\n=== Chapter 3 complete ===")