# ============================================================
#  HOOFDSTUK 1: Kinetische validatie Model 1
#  Vier analyses:
#   1a. Één-stap groeicurve: latente periode en burst size
#   1c. Gevoeligheid voor alfa_ads en tau
#   1d. One-step growth experiment: 5 MOI-curves + biomassa
#   1e. Groei op glucose vs maltose vergelijking
#  (1b verwijderd — MOI-sweep niet langer opgenomen)
# ============================================================
include("./Bin_model_1/dFBA.jl")
include("./Bin_model_1/parameters.jl")
include("./Bin_model_1/FBA.jl")

using Plots, Statistics, SciMLBase
using COBREXA, AbstractFBCModels
using DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

model_path = joinpath(@__DIR__, "iJO1366.xml")
model      = loadFBAmodel(model_path)

alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 79/60;  b = 170.0   # tau consistent met TAU_L = 79/60
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [7.24, 2.26, 0.0, 10.0]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
h_release      = 1.71e-12;  duration = 20.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
alfa_ads       = 1.3e-9
fbaCache       = buildFbaCache(model, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
E_coli_cellDW  = 1.0e-12

function make_p1(N0, t_inf, moi; tau_=tau, alfa_=alfa_ads, vmax_=V_max)
    Parameters(duration, N0, alpha_syn, beta_deg, K_s, vmax_, p_pref,
        tau_, b, alfa_, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        fbaCache, 0.0, zeros(4), t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0, mu_max_vector, e_max_vector, 0.0, 0.0)
end

# ============================================================
#  1a. Growth curve validation (lage MOI = lytisch regime)
# ============================================================
println("=== 1a: Growth curve validation ===")
p_ref = make_p1(1e9, 2.0, 0.001)
sol   = run(p_ref)

t   = sol.t
S   = [sol.u[i][Sind_S] for i in eachindex(sol.u)]
I   = [sol.u[i][Sind_I] for i in eachindex(sol.u)]
L   = [sol.u[i][Sind_L] for i in eachindex(sol.u)]
P   = [sol.u[i][Pfind]  for i in eachindex(sol.u)]
MOI_raw = [begin
    s     = sol.u[i][Sind_S]
    p_val = sol.u[i][Pfind]
    (s > CELL_THRESHOLD && p_val > PHAGE_THRESHOLD) ?
        min(p_val / s, 10.0) : NaN
end for i in eachindex(sol.u)]

MOI = copy(MOI_raw)
last_valid = isnan(MOI_raw[1]) ? 0.0 : MOI_raw[1]
for i in eachindex(MOI_raw)
    if !isnan(MOI_raw[i])
        last_valid = MOI_raw[i]
    else
        MOI[i] = last_valid
    end
end
X   = S .+ I .+ L

i_max_I = argmax(I)
t_lysis = t[i_max_I]
println("  Lysis peak at t = $(round(t_lysis,digits=2)) h (expected ≈ $(2.0 + 1.32) h)")

# Verwijder size uit de individuele subplots pa, pb, pc_moi
# en zet alleen size op de gecombineerde figuur

# Fix: gebruik gr() expliciet en zet dpi
gr()

S_plot  = max.(S, 1.0)
I_plot  = max.(I, 1.0)
L_plot  = max.(L, 1.0)
X_plot  = max.(X, 1.0)
P_plot  = max.(P, 1.0)
MOI_plot = max.(MOI .+ 1e-6, 1e-6)

pa = plot(t, S_plot,  label="Susceptible S",  color=:blue,  lw=2)
plot!(pa, t, I_plot,  label="Lytic I",        color=:red,   lw=2)
plot!(pa, t, L_plot,  label="Lysogenic L",    color=:green, lw=2)
plot!(pa, t, X_plot,  label="Total X",        color=:black, lw=2)
plot!(pa,
    yscale=:log10,
    ylims=(10.0, 1e12),
    xlabel="t [h]", ylabel="Cells L⁻¹",
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    legend=:bottomleft,
    bottom_margin=5Plots.mm, left_margin=8Plots.mm)

pb = plot(t, P_plot, label="Free phages P", color=:darkred, lw=2)
plot!(pb,
    yscale=:log10,
    ylims=(10.0, 1e13),
    xlabel="t [h]", ylabel="Phages L⁻¹",
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    legend=:topleft,
    bottom_margin=5Plots.mm, left_margin=8Plots.mm)

pc_moi = plot(t, MOI_plot, label="MOI", color=:purple, lw=2)
plot!(pc_moi,
    yscale=:log10,
    ylims=(1e-6, 20.0),
    xlabel="t [h]", ylabel="MOI [-]",
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    legend=:topleft,
    bottom_margin=5Plots.mm, left_margin=8Plots.mm)

fig1a = plot(pa, pb, pc_moi,
    layout=(1,3),
    size=(1500, 500),
    left_margin=2Plots.mm,
    bottom_margin=2Plots.mm,
    top_margin=2Plots.mm,
    right_margin=2Plots.mm)
savefig(fig1a, "h1a_groeicurve_validatie.png")
println("  Figure saved: h1a_groeicurve_validatie.png")


# ============================================================
#  1c. Sensitivity analysis: alfa_ads and tau
# ============================================================
println("\n=== 1c: Sensitivity analysis ===")

alfa_waarden = [1e-11, 1e-10, 5e-10, 1e-9]
# tau sensitivity: varieer TAU_L via p.tau, TAU_D schaalt mee als 28/79 * tau_
tau_waarden  = [0.5, 1.0, 79/60, 2.0]

max_P_alfa = Float64[]
for alfa in alfa_waarden
    sol_a = run(make_p1(1e9, 2.0, 0.001; alfa_=alfa))
    push!(max_P_alfa, SciMLBase.successful_retcode(sol_a) ?
        maximum([sol_a.u[i][Pfind] for i in eachindex(sol_a.u)]) : 0.0)
    println("  alfa_ads=$alfa: max_P=$(round(max_P_alfa[end], sigdigits=2))")
end

t_piek_tau = Float64[]
for tau_ in tau_waarden
    sol_t = run(make_p1(1e9, 2.0, 0.001; tau_=tau_))
    if SciMLBase.successful_retcode(sol_t)
        I_t = [sol_t.u[i][Sind_I] for i in eachindex(sol_t.u)]
        t_t = sol_t.t
        idx = argmax(I_t)
        println("  tau=$(round(tau_, digits=3)) h:")
        println("    max I     = $(maximum(I_t))")
        println("    argmax I  = $idx")
        println("    t_piek    = $(t_t[idx])")
        push!(t_piek_tau, t_t[idx])
    else
        println("  tau=$(round(tau_, digits=3)) h: SOLVER FAILED")
        push!(t_piek_tau, NaN)
    end
end

pf = plot(alfa_waarden, max_P_alfa,
    xlabel="α_ads [L cell⁻¹ h⁻¹]", ylabel="Max phages L⁻¹",
    xscale=:log10,
    color=:steelblue, lw=2, marker=:circle, markersize=7, legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=8Plots.mm,
    size=(500, 400))
pg = plot(tau_waarden, t_piek_tau,
    xlabel="τ [h]", ylabel="t_peak [h]",
    color=:darkorange, lw=2, marker=:circle, markersize=7, legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=8Plots.mm,
    size=(500, 400))
fig1c = plot(pf, pg, layout=(1,2), size=(1000,420), margin=6Plots.mm)
savefig(fig1c, "h1c_gevoeligheidsanalyse.png")
println("  Figure saved: h1c_gevoeligheidsanalyse.png")

# ============================================================
#  1d. One-step growth experiment: 5 MOI-waarden op één grafiek
# ============================================================
println("\n=== 1d: One-step growth experiment (5 MOI curves) ===")

moi_curves   = [0.001, 0.01, 0.1, 1.0, 5.0]
kleuren      = [:darkred :red :orange :steelblue :darkblue]
labels_P     = ["MOI=0.001" "MOI=0.01" "MOI=0.1" "MOI=1.0" "MOI=5.0"]

fig1d_P  = plot(xlabel="t [h]", ylabel="Free phages L⁻¹",
    yscale=:log10, ylims=(1, :auto), legend=:topleft, size=(650,420),
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
fig1d_X  = plot(xlabel="t [h]", ylabel="Total biomass (cells L⁻¹)",
    yscale=:log10, ylims=(1e4, :auto), legend=:bottomleft, size=(650,420),
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)

for (idx, moi) in enumerate(moi_curves)
    sol_d = run(make_p1(1e9, 2.0, moi))
    if SciMLBase.successful_retcode(sol_d)
        t_d = sol_d.t
        P_d = max.([sol_d.u[i][Pfind]  for i in eachindex(sol_d.u)], 1.0)
        S_d = [sol_d.u[i][Sind_S] for i in eachindex(sol_d.u)]
        I_d = [sol_d.u[i][Sind_I] for i in eachindex(sol_d.u)]
        L_d = [sol_d.u[i][Sind_L] for i in eachindex(sol_d.u)]
        X_d = max.(S_d .+ I_d .+ L_d, 1.0)
        plot!(fig1d_P, t_d, P_d, label=labels_P[idx], color=kleuren[idx], lw=2)
        plot!(fig1d_X, t_d, X_d, label=labels_P[idx], color=kleuren[idx], lw=2)
    end
end

vline!(fig1d_P, [2.0],        color=:gray,  lw=1, linestyle=:dash, label="t_inf")
vline!(fig1d_P, [2.0 + 1.32], color=:black, lw=1, linestyle=:dot,  label="t_inf + τ")
vline!(fig1d_X, [2.0],        color=:gray,  lw=1, linestyle=:dash, label="t_inf")

fig1d = plot(fig1d_P, fig1d_X, layout=(1,2), size=(1100,430), margin=6Plots.mm)
savefig(fig1d, "h1d_onestep_growth.png")
println("  Figure saved: h1d_onestep_growth.png")

# ============================================================
#  1e. Growth on glucose vs maltose
# ============================================================
println("\n=== 1e: Growth on glucose vs maltose ===")

V_max_glc  = [12.7, 0.0, 0.0, 4.0]
V_max_malt = [0.0, 3.75, 0.0, 4.0]
V_max_both = [12.7, 3.75, 0.0, 4.0]

function run_groei(vmax_scenario, glc0, malt0)
    p = make_p1(1e9, 999.0, 0.0; vmax_=vmax_scenario)
    u0 = zeros(13)
    u0[Sind]   = [glc0, malt0, 5.42, 0.0]
    u0[Eind]   = [0.95, 0.01, 0.01, 0.01]
    u0[Sind_S] = 1e9
    tspan = (0.0, 10.0)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:10.0)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
    function fbaAffect!(integrator)
        fbaUpdate!(integrator.u, p)
        enforceThreshold!(integrator.u)
    end
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator)  = any(x -> x < 0.0, u)
    domainAffect!(integrator)          = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    prob = DDEProblem(simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(prob, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=fbaUpdateTimepoints,
        callback=CallbackSet(domainCallBack, infectionCallBack, fbaCallBack))
end

sol_glc  = run_groei(V_max_glc,  4.44, 0.0)
sol_malt = run_groei(V_max_malt, 0.0,  2.337)
sol_both = run_groei(V_max_both, 4.44, 2.337)

function extract_groei(sol)
    t_g  = sol.t
    N_g  = [sol.u[i][Sind_S] for i in eachindex(sol.u)]
    glc  = [sol.u[i][Sind[1]] for i in eachindex(sol.u)]
    malt = [sol.u[i][Sind[2]] for i in eachindex(sol.u)]
    e_mal= [sol.u[i][Eind[2]] for i in eachindex(sol.u)]
    return t_g, N_g, glc, malt, e_mal
end

t_g, N_g, glc_g, malt_g, emal_g = extract_groei(sol_glc)
t_m, N_m, glc_m, malt_m, emal_m = extract_groei(sol_malt)
t_b, N_b, glc_b, malt_b, emal_b = extract_groei(sol_both)

function schat_mu(t, N)
    idx_start = findfirst(N .> 1e6)
    idx_end   = isnothing(idx_start) ? length(N) : min(idx_start+50, length(N))
    if !isnothing(idx_start) && N[idx_end] > N[idx_start] && t[idx_end] > t[idx_start]
        return log(N[idx_end]/N[idx_start]) / (t[idx_end]-t[idx_start])
    end
    return NaN
end

mu_glc  = schat_mu(t_g, N_g)
mu_malt = schat_mu(t_m, N_m)
mu_both = schat_mu(t_b, N_b)

println("  Estimated µ glucose-only  : $(round(mu_glc,  digits=3)) h⁻¹")
println("  Estimated µ maltose-only  : $(round(mu_malt, digits=3)) h⁻¹")
println("  Estimated µ glucose+malt  : $(round(mu_both, digits=3)) h⁻¹")

p1e_N = plot(xlabel="t [h]", ylabel="Number of cells L⁻¹",
    yscale=:log10, ylims=(1e7, :auto), legend=:topleft,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p1e_N, t_g, max.(N_g, 1.0), label="Glucose only",      color=:steelblue,  lw=2.5)
plot!(p1e_N, t_m, max.(N_m, 1.0), label="Maltose only",      color=:darkorange, lw=2.5, linestyle=:dash)
plot!(p1e_N, t_b, max.(N_b, 1.0), label="Glucose + Maltose", color=:darkgreen,  lw=2.5, linestyle=:dot)
annotate!(p1e_N, [(8.5, maximum(filter(x->x>1, N_g))*0.6,
    text("µ≈$(round(mu_glc,digits=2)) h⁻¹", 9, :steelblue, :left))])
annotate!(p1e_N, [(8.5, maximum(filter(x->x>1, N_m))*0.6,
    text("µ≈$(round(mu_malt,digits=2)) h⁻¹", 9, :darkorange, :left))])
annotate!(p1e_N, [(8.5, maximum(filter(x->x>1, N_b))*1.1,
    text("µ≈$(round(mu_both,digits=2)) h⁻¹", 9, :darkgreen, :left))])

p1e_S = plot(xlabel="t [h]", ylabel="Concentration [mmol L⁻¹]",
    legend=:topright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p1e_S, t_b, glc_b,  label="Glucose",  color=:steelblue,  lw=2)
plot!(p1e_S, t_b, malt_b, label="Maltose",  color=:darkorange, lw=2, linestyle=:dash)

p1e_emal = plot(xlabel="t [h]", ylabel="e_mal [-]",
    legend=:topleft,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p1e_emal, t_g, emal_g, label="Glucose only",      color=:steelblue,  lw=2.5)
plot!(p1e_emal, t_m, emal_m, label="Maltose only",      color=:darkorange, lw=2.5, linestyle=:dash)
plot!(p1e_emal, t_b, emal_b, label="Glucose + Maltose", color=:darkgreen,  lw=2.5, linestyle=:dot)

fig1e = plot(p1e_N, p1e_S, p1e_emal, layout=(1,3), size=(1300,430), margin=6Plots.mm)
savefig(fig1e, "h1e_glucose_vs_maltose.png")
println("  Figure saved: h1e_glucose_vs_maltose.png")

println("\n=== Chapter 1 complete ===")