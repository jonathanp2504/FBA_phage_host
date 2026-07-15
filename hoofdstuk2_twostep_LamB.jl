# ============================================================
#  CHAPTER 4: One-step and two-step phage adsorption dynamics
#
#  No production or toxicity anywhere in this script (f_prod = 0,
#  Y_benz = 0, k_tox = 0, beta_benz = 0 throughout): benzonase is
#  only introduced in Chapter 5.
#
#  IMPORTANT: the one-step model here uses the original
#  Bin_model_1/dFBA.jl + parameters.jl + FBA.jl implementation,
#  NOT the later Model1.jl module. Model1.jl collapses the decision
#  delay tau_d and the latent period tau_l into a single delay
#  (h(p, t - p.tau) used for both the decision and the lysis
#  history lookup), whereas the original implementation correctly
#  uses two separate delays (tau_d = 28/60 h fixed, tau_l = p.tau).
#  That collapse was producing a much larger, unrealistic biomass
#  drop upon infection in the one-step simulations. Model2.jl is
#  unaffected (it already uses two separate delayed lookups) and is
#  kept as-is for the two-step figures.
#
#  Figures produced (matching Chapter 4's figure references):
#   h1a_groeicurve_validatie.png    - one-step population dynamics, MOI=0.001
#   h1d_onestep_growth.png          - one-step growth experiment, 5 MOI values
#   h1_onestep_crash.png            - one-step at the "biologically correct" k_ads
#   h2_twostep_highrate_works.png  - two-step at the same rate, stable
#   h2a_model2_moi_vergelijking.png - two-step population dynamics, 4 MOI values
#   h2a_model2_moi_fagen.png        - two-step phage pools, 4 MOI values
#   h2b_onestep_vs_twostep_alfa.png - appendix: one-step vs two-step comparison
# ============================================================
include("./Bin_model_1/dFBA.jl")
include("./Bin_model_1/parameters.jl")
include("./Bin_model_1/FBA.jl")
include("Model2.jl")
using Plots, Statistics, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

gr()

# ------------------------------------------------------------
#  Shared model and metabolic parameters
# ------------------------------------------------------------
model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;   beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 1.32;  b = 170.0                 # latent period [h], burst size
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max_malt     = [0.0, 2.26, 0.0, 10.0]            # maltose-only growth (two-step figures)
V_max_full     = [7.24, 2.26, 0.0, 10.0]           # glucose+maltose growth (one-step figures)
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 1.71e-12
duration       = 20.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW  = 1.0e-12

k_ads_low   = 1.3e-9    # one-step reference adsorption constant
k_att_high  = 7.92e-8   # two-step attachment constant (also used as the "crash" k_ads)

model_one_step = loadFBAmodel(model_path)
fbaCache1      = buildFbaCache(model_one_step, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")

baseModel2  = Model2.loadFBAmodel(model_path)
all_ex_ids  = [id for id in keys(baseModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2   = Model2.buildFbaCache(baseModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2 = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                    "R_BIOMASS_Ec_iJO1366_core_53p95M")

# ------------------------------------------------------------
#  Parameter constructors -- no production/toxicity (all zero)
# ------------------------------------------------------------
function make_p1(N0, t_inf, moi; k_ads=k_ads_low, vmax=V_max_full, tau_=tau)
    Parameters(duration, N0, alpha_syn, beta_deg, K_s,
        vmax, p_pref, tau_, b, k_ads, E_coli_cellDW,
        MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        fbaCache1, 0.0, zeros(4), t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0, mu_max_vector, e_max_vector, 0.0, 0.0)
end

function make_p2(N0, t_inf, moi; k_att=k_att_high, vmax=V_max_malt)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, vmax, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        k_att, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0,
        mu_max_vector, e_max_vector, 0.0, 0.0)
end

# ------------------------------------------------------------
#  Shared, minimal plot style (short labels, legend outside the
#  axes so it never overlaps a curve, generous margins)
# ------------------------------------------------------------
default_style = (tickfontsize=10, guidefontsize=11, legendfontsize=8,
                  bottom_margin=8Plots.mm, left_margin=10Plots.mm, top_margin=3Plots.mm)

# ============================================================
#  h1a: one-step population dynamics at MOI = 0.001
# ============================================================
println("=== h1a: one-step population dynamics (MOI=0.001) ===")
sol_h1a = run(make_p1(1e9, 2.0, 0.001))

t_a = sol_h1a.t
S_a = [sol_h1a.u[i][Sind_S] for i in eachindex(sol_h1a.u)]
I_a = [sol_h1a.u[i][Sind_I] for i in eachindex(sol_h1a.u)]
L_a = [sol_h1a.u[i][Sind_L] for i in eachindex(sol_h1a.u)]
P_a = [sol_h1a.u[i][Pfind]  for i in eachindex(sol_h1a.u)]
X_a = S_a .+ I_a .+ L_a

MOI_raw = [begin
    s = sol_h1a.u[i][Sind_S]; p_ = sol_h1a.u[i][Pfind]
    (s > CELL_THRESHOLD && p_ > PHAGE_THRESHOLD) ? min(p_/s, 10.0) : NaN
end for i in eachindex(sol_h1a.u)]
MOI_a = copy(MOI_raw)
last_valid = isnan(MOI_raw[1]) ? 0.0 : MOI_raw[1]
for i in eachindex(MOI_raw)
    MOI_a[i] = isnan(MOI_raw[i]) ? last_valid : MOI_raw[i]
    last_valid = MOI_a[i]
end

p_h1a_cells = plot(t_a, [max.(S_a,1.0) max.(I_a,1.0) max.(L_a,1.0) max.(X_a,1.0)],
    label=["Naive" "Lytic" "Lysogenic" "Total"],
    color=[:blue :red :green :black], lw=2,
    yscale=:log10, ylims=(10.0,1e12),
    xlabel="t [h]", ylabel="Cells/L",
    legend=:outertop, legend_column=-1; default_style...)

p_h1a_phage = plot(t_a, max.(P_a,1.0), label="Free phages", color=:darkred, lw=2,
    yscale=:log10, ylims=(10.0,1e13),
    xlabel="t [h]", ylabel="Phages/L", legend=false; default_style...)

p_h1a_moi = plot(t_a, max.(MOI_a,1e-6), label="MOI", color=:purple, lw=2,
    yscale=:log10, ylims=(1e-6,20.0),
    xlabel="t [h]", ylabel="MOI [-]", legend=false; default_style...)

fig_h1a = plot(p_h1a_cells, p_h1a_phage, p_h1a_moi, layout=(1,3), size=(1650,480), margin=4Plots.mm)
savefig(fig_h1a, "h1a_groeicurve_validatie.png")
println("  Figure saved: h1a_groeicurve_validatie.png")

# ============================================================
#  h1d: one-step growth experiment, 5 MOI values
# ============================================================
println("\n=== h1d: one-step growth experiment (5 MOI values) ===")
moi_curves = [0.001, 0.01, 0.1, 1.0, 5.0]
colors_d   = [:darkred, :red, :orange, :steelblue, :darkblue]

p_h1d_phage = plot(xlabel="t [h]", ylabel="Phages/L",
    yscale=:log10, ylims=(1,:auto),
    legend=:outertop, legend_column=-1; default_style...)
p_h1d_bio = plot(xlabel="t [h]", ylabel="Biomass [cells/L]",
    yscale=:log10, ylims=(1e4,:auto), legend=false; default_style...)

for (idx, moi) in enumerate(moi_curves)
    sol_d = run(make_p1(1e9, 2.0, moi))
    if SciMLBase.successful_retcode(sol_d)
        t_d = sol_d.t
        P_d = max.([sol_d.u[i][Pfind] for i in eachindex(sol_d.u)], 1.0)
        S_d = [sol_d.u[i][Sind_S] for i in eachindex(sol_d.u)]
        I_d = [sol_d.u[i][Sind_I] for i in eachindex(sol_d.u)]
        L_d = [sol_d.u[i][Sind_L] for i in eachindex(sol_d.u)]
        X_d = max.(S_d .+ I_d .+ L_d, 1.0)
        plot!(p_h1d_phage, t_d, P_d, label="MOI=$moi", color=colors_d[idx], lw=2)
        plot!(p_h1d_bio,   t_d, X_d, label="MOI=$moi", color=colors_d[idx], lw=2)
    end
end
vline!(p_h1d_phage, [2.0, 3.32], color=[:gray :black], lw=1, linestyle=[:dash :dot], label=["t_inf" "t_inf+tau_l"])
vline!(p_h1d_bio,   [2.0], color=:gray, lw=1, linestyle=:dash, label="")

fig_h1d = plot(p_h1d_phage, p_h1d_bio, layout=(1,2), size=(1200,460), margin=5Plots.mm)
savefig(fig_h1d, "h1d_onestep_growth.png")
println("  Figure saved: h1d_onestep_growth.png")

# ============================================================
#  h1_onestep_crash: one-step model at the "biologically correct"
#  (higher) adsorption rate -- expected to become numerically unstable
# ============================================================
println("\n=== h1_onestep_crash: one-step at k_ads = k_att ===")
sol_crash = run(make_p1(1e9, 2.0, 0.01; k_ads=k_att_high))
t_c = sol_crash.t
S_c = [sol_crash.u[i][Sind_S] for i in eachindex(sol_crash.u)]
I_c = [sol_crash.u[i][Sind_I] for i in eachindex(sol_crash.u)]
L_c = [sol_crash.u[i][Sind_L] for i in eachindex(sol_crash.u)]
P_c = [sol_crash.u[i][Pfind]  for i in eachindex(sol_crash.u)]
X_c = max.(S_c .+ I_c .+ L_c, 1.0)

p_crash_bio = plot(t_c, X_c, label="Total biomass", color=:black, lw=2,
    yscale=:log10, ylims=(1,:auto),
    xlabel="t [h]", ylabel="Biomass [cells/L]", legend=false; default_style...)
p_crash_phage = plot(t_c, max.(P_c,1.0), label="Free phages", color=:darkred, lw=2,
    yscale=:log10, ylims=(1,:auto),
    xlabel="t [h]", ylabel="Phages/L", legend=false; default_style...)

fig_crash = plot(p_crash_bio, p_crash_phage, layout=(1,2), size=(1200,460), margin=5Plots.mm)
savefig(fig_crash, "h1_onestep_crash.png")
println("  Figure saved: h1_onestep_crash.png")
println("  Solver return code: $(sol_crash.retcode)")

# ============================================================
#  h2_twostep_highrate_works: two-step model at the same high rate,
#  shown to remain stable
# ============================================================
println("\n=== h2_twostep_highrate_works: two-step at k_att = $(k_att_high) ===")
sol_works = Model2.run(make_p2(1e9, 2.0, 0.01))
t_w = sol_works.t
N_w = [sol_works.u[i][Model2.Nind] for i in eachindex(sol_works.u)]
Lyt_w = [sol_works.u[i][Model2.Lind] for i in eachindex(sol_works.u)]
Lys_w = [sol_works.u[i][Model2.lind] for i in eachindex(sol_works.u)]
P_w = [sol_works.u[i][Model2.Pfind] for i in eachindex(sol_works.u)]
X_w = max.(N_w .+ Lyt_w .+ Lys_w, 1.0)

p_works_bio = plot(t_w, X_w, label="Total biomass", color=:black, lw=2,
    yscale=:log10, ylims=(1e4,:auto),
    xlabel="t [h]", ylabel="Biomass [cells/L]", legend=false; default_style...)
p_works_phage = plot(t_w, max.(P_w,1.0), label="Free phages", color=:darkred, lw=2,
    yscale=:log10, ylims=(1,:auto),
    xlabel="t [h]", ylabel="Phages/L", legend=false; default_style...)

fig_works = plot(p_works_bio, p_works_phage, layout=(1,2), size=(1200,460), margin=5Plots.mm)
savefig(fig_works, "h2_twostep_highrate_works.png")
println("  Figure saved: h2_twostep_highrate_works.png")
println("  Solver return code: $(sol_works.retcode)")

# ============================================================
#  h2a: two-step population dynamics at 4 MOI values
# ============================================================
println("\n=== h2a: two-step population dynamics (4 MOI values) ===")
moi_compare = [0.001, 0.1, 1.0, 5.0]
panels_cel  = []
panels_faag = []

for (idx, moi) in enumerate(moi_compare)
    sol2 = Model2.run(make_p2(1e9, 2.0, moi))
    t2   = sol2.t
    N2   = [sol2.u[i][Model2.Nind]  for i in eachindex(sol2.u)]
    Lyt2 = [sol2.u[i][Model2.Lind]  for i in eachindex(sol2.u)]
    Lys2 = [sol2.u[i][Model2.lind]  for i in eachindex(sol2.u)]
    Pf2  = max.([sol2.u[i][Model2.Pfind] for i in eachindex(sol2.u)], 1.0)
    Pa2  = max.([sol2.u[i][Model2.Paind] for i in eachindex(sol2.u)], 1.0)

    p_cel = plot(t2, [max.(N2,1.0) max.(Lyt2,1.0) max.(Lys2,1.0)],
        label = idx==1 ? ["Naive" "Lytic" "Lysogenic"] : false,
        color=[:blue :red :green], lw=2, yscale=:log10, ylims=(1,:auto),
        title="MOI = $moi",
        xlabel="t [h]", ylabel="Cells/L",
        legend = idx==1 ? :outertop : false, legend_column=-1; default_style...)
    push!(panels_cel, p_cel)

    p_faag = plot(t2, Pf2, label = idx==1 ? "Free" : false, color=:darkred, lw=2,
        yscale=:log10, ylims=(1,:auto), title="MOI = $moi",
        xlabel="t [h]", ylabel="Phages/L",
        legend = idx==1 ? :outertop : false, legend_column=-1; default_style...)
    plot!(p_faag, t2, Pa2, label = idx==1 ? "Attached" : false,
        color=:purple, lw=1.5, linestyle=:dash)
    push!(panels_faag, p_faag)
end

fig2a_cel  = plot(panels_cel...,  layout=(2,2), size=(1100,760), margin=4Plots.mm)
fig2a_faag = plot(panels_faag..., layout=(2,2), size=(1100,760), margin=4Plots.mm)
savefig(fig2a_cel,  "h2a_model2_moi_vergelijking.png")
savefig(fig2a_faag, "h2a_model2_moi_fagen.png")
println("  Figure saved: h2a_model2_moi_vergelijking.png")
println("  Figure saved: h2a_model2_moi_fagen.png")

# ============================================================
#  h2b (appendix): one-step (low k_ads), one-step (k_ads = k_att),
#  and two-step model compared at identical parameters
# ============================================================
println("\n=== h2b: one-step vs two-step comparison (appendix) ===")
moi_test = 0.01

sol_m1_low  = run(make_p1(1e9, 2.0, moi_test; k_ads=k_ads_low))
sol_m1_high = run(make_p1(1e9, 2.0, moi_test; k_ads=k_att_high))
sol_m2_ref  = Model2.run(make_p2(1e9, 2.0, moi_test))

function extract_m1(sol)
    t = sol.t
    S = [sol.u[i][Sind_S] for i in eachindex(sol.u)]
    I = [sol.u[i][Sind_I] for i in eachindex(sol.u)]
    L = [sol.u[i][Sind_L] for i in eachindex(sol.u)]
    P = max.([sol.u[i][Pfind] for i in eachindex(sol.u)], 1.0)
    X = S .+ I .+ L
    return (t=t, I=I, L=L, P=P, X=X)
end
function extract_m2(sol)
    t   = sol.t
    N   = [sol.u[i][Model2.Nind]  for i in eachindex(sol.u)]
    Lyt = [sol.u[i][Model2.Lind]  for i in eachindex(sol.u)]
    Lys = [sol.u[i][Model2.lind]  for i in eachindex(sol.u)]
    P   = max.([sol.u[i][Model2.Pfind] for i in eachindex(sol.u)], 1.0)
    X   = N .+ Lyt .+ Lys
    return (t=t, Lyt=Lyt, Lys=Lys, P=P, X=X)
end

r1_low  = extract_m1(sol_m1_low)
r1_high = extract_m1(sol_m1_high)
r2_ref  = extract_m2(sol_m2_ref)

lys_frac(I, L) = [ (I[i]+L[i]) > 1.0 ? L[i]/(I[i]+L[i]) : 0.0 for i in eachindex(I) ]

p2b_X = plot(xlabel="t [h]", ylabel="Biomass [cells/L]",
    yscale=:log10, ylims=(1e4,:auto),
    legend=:outertop, legend_column=-1; default_style...)
plot!(p2b_X, r1_low.t,  max.(r1_low.X,1.0),  label="One-step (low)",  color=:steelblue, lw=2)
plot!(p2b_X, r1_high.t, max.(r1_high.X,1.0), label="One-step (high)", color=:red, lw=2)
plot!(p2b_X, r2_ref.t,  max.(r2_ref.X,1.0),  label="Two-step",        color=:darkgreen, lw=2)

p2b_P = plot(xlabel="t [h]", ylabel="Phages/L",
    yscale=:log10, ylims=(1,:auto),
    legend=:outertop, legend_column=-1; default_style...)
plot!(p2b_P, r1_low.t,  r1_low.P,  label="One-step (low)",  color=:steelblue, lw=2)
plot!(p2b_P, r1_high.t, r1_high.P, label="One-step (high)", color=:red, lw=2)
plot!(p2b_P, r2_ref.t,  r2_ref.P,  label="Two-step",        color=:darkgreen, lw=2)

p2b_frac = plot(xlabel="t [h]", ylabel="Lysogenic fraction [-]",
    ylims=(0,1.05),
    legend=:outertop, legend_column=-1; default_style...)
plot!(p2b_frac, r1_low.t,  lys_frac(r1_low.I,  r1_low.L),  label="One-step (low)",  color=:steelblue, lw=2)
plot!(p2b_frac, r1_high.t, lys_frac(r1_high.I, r1_high.L), label="One-step (high)", color=:red, lw=2)
plot!(p2b_frac, r2_ref.t,  lys_frac(r2_ref.Lyt, r2_ref.Lys), label="Two-step",      color=:darkgreen, lw=2)

fig2b = plot(p2b_X, p2b_P, p2b_frac, layout=(1,3), size=(1650,480), margin=4Plots.mm)
savefig(fig2b, "h2b_onestep_vs_twostep_alfa.png")
println("  Figure saved: h2b_onestep_vs_twostep_alfa.png")

println("\n=== Chapter 4 simulations complete ===")