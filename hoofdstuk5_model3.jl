# ============================================================
#  HOOFDSTUK 5 — Model 3: Process optimisation
# ============================================================
include("Model3.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

const BENZ_REF_3b = 0.124

model_path    = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn     = 2.0;  beta_deg = 0.5
K_s           = [0.0061, 9.4e-4, 0.0543, 8.33]
tau           = 1.32
p_pref        = [0.8925, 0.08925, 0.008925, 0.008925]
V_max         = [0.0, 2.26, 0.0, 10.0]
exchange_ids  = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids= ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                 "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                 "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                 "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                 "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values     = [180.16, 342.3, 92.09, 60.05]
h_release     = 1.71e-12;  duration = 40.0
mu_max_vector = [0.76, 0.76, 1.10, 0.30]
e_max_vector  = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW = 1.0e-12

naiveModel3   = Model3.loadFBAmodel(model_path)
lysogenModel3 = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
all_ex_ids    = [id for id in keys(naiveModel3.reactions) if startswith(id, "R_EX_")]
naiveFba3     = Model3.buildFbaCache(naiveModel3,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3   = Model3.buildFbaCache(lysogenModel3, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")

function make_p3b(N0, t_inf, moi, k_tox=0.35, f_prod=0.0015)
    Model3.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, 170.0, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf, moi*N0,
        "R_BENZ_prod", k_tox, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod)
end

# ============================================================
#  5a. Optimal run Model 3
# ============================================================
println("=== 5a (M3): Optimal run Model 3 ===")

opt_moi      = 0.0094; opt_t_inf = 3.77; opt_biomassa = 1e9

p_opt   = make_p3b(opt_biomassa, opt_t_inf, opt_moi)
sol_opt = Model3.run(p_opt)
t_opt   = sol_opt.t
B_opt   = [sol_opt.u[i][Model3.Benzind] for i in eachindex(sol_opt.u)]
Pf_opt  = [sol_opt.u[i][Model3.Pfind]  for i in eachindex(sol_opt.u)]
N_opt   = [sol_opt.u[i][Model3.Nind]   for i in eachindex(sol_opt.u)]
PN_opt  = [N_opt[i]>1.0 ? Pf_opt[i]/N_opt[i] : 0.0 for i in eachindex(t_opt)]

pa = plot(t_opt, B_opt, label="Benzonase", color=:green, lw=2, fill=(0,0.2,:green),
    xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]", legend=:topleft,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
pb = plot(t_opt, max.(PN_opt, 1e-10), label="P/N", color=:red, lw=2, yscale=:log10,
    xlabel="t [h]", ylabel="P/N [-]", legend=:topright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
fig5a_m3 = plot(pa, pb, layout=(1,2), size=(900,430), margin=6Plots.mm)
savefig(fig5a_m3, "h5a_m3_optimalisatie.png")
println("  Saved: h5a_m3_optimalisatie.png")

# ============================================================
#  5_pareto. Pareto punten A, B, C — gedetailleerde tijdseries
#
#  A: MOI=3.0,  t_inf=0.5h  → maximale benzonase, lange dt
#  C: MOI=0.01, t_inf=3.0h  → compromispunt (knie van het front)
#  B: MOI=0.1,  t_inf=7.0h  → kortste productietijd, lage benzonase
# ============================================================
println("\n=== 5_pareto: Pareto punten A, B, C tijdseries ===")

pareto_punten = [
    (label="A", moi=3.0,  t_inf=0.5, kleur=:darkred),
    (label="C", moi=0.01, t_inf=3.0, kleur=:steelblue),
    (label="B", moi=0.1,  t_inf=7.0, kleur=:darkgreen),
]

panels_benz = []; panels_pop = []; panels_pf = []

for pt in pareto_punten
    sol_pt = Model3.run(make_p3b(1e9, pt.t_inf, pt.moi))

    t_pt  = sol_pt.t
    B_pt  = [sol_pt.u[i][Model3.Benzind] for i in eachindex(sol_pt.u)]
    N_pt  = [sol_pt.u[i][Model3.Nind]    for i in eachindex(sol_pt.u)]
    D_pt  = [sol_pt.u[i][Model3.Dind]    for i in eachindex(sol_pt.u)]
    L_pt  = [sol_pt.u[i][Model3.Lind]    for i in eachindex(sol_pt.u)]
    l_pt  = [sol_pt.u[i][Model3.lind]    for i in eachindex(sol_pt.u)]
    Pf_pt = [max(sol_pt.u[i][Model3.Pfind], 1.0) for i in eachindex(sol_pt.u)]

    # Benzonase
    p_b = plot(t_pt, B_pt,
        label="Benzonase", color=pt.kleur, lw=2,
        fill=(0, 0.15, pt.kleur),
        xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
        title="Point $(pt.label) — MOI=$(pt.moi), t_inf=$(pt.t_inf)h",
        titlefontsize=11,
        tickfontsize=10, guidefontsize=12, legendfontsize=9,
        legend=:topleft,
        bottom_margin=6Plots.mm, left_margin=10Plots.mm)
    vline!(p_b, [pt.t_inf], color=:gray, lw=1, linestyle=:dash, label="t_inf")
    push!(panels_benz, p_b)

    # Populaties
    p_pop = plot(t_pt, max.(N_pt, 1.0),
        label="Naive N", color=:blue, lw=2,
        yscale=:log10, ylims=(1, :auto),
        xlabel="t [h]", ylabel="Cells L⁻¹",
        title="Point $(pt.label) — populations",
        titlefontsize=11,
        tickfontsize=10, guidefontsize=12, legendfontsize=9,
        legend=:outertopright,
        bottom_margin=6Plots.mm, left_margin=10Plots.mm)
    plot!(p_pop, t_pt, max.(L_pt, 1.0), label="Lytic L",   color=:red,   lw=2)
    plot!(p_pop, t_pt, max.(l_pt, 1.0), label="Lysogen l", color=:green, lw=2)
    vline!(p_pop, [pt.t_inf], color=:gray, lw=1, linestyle=:dash, label="t_inf")
    push!(panels_pop, p_pop)

    # Fagen
    p_pf = plot(t_pt, Pf_pt,
        label="Free phages", color=:darkred, lw=2,
        yscale=:log10, ylims=(1, :auto),
        xlabel="t [h]", ylabel="Phages L⁻¹",
        title="Point $(pt.label) — phages",
        titlefontsize=11,
        tickfontsize=10, guidefontsize=12, legendfontsize=9,
        legend=:topleft,
        bottom_margin=6Plots.mm, left_margin=10Plots.mm)
    vline!(p_pf, [pt.t_inf], color=:gray, lw=1, linestyle=:dash, label="t_inf")
    push!(panels_pf, p_pf)
end

# Drie rijen (A, C, B) × drie kolommen (benzonase, populaties, fagen)
fig5_pareto = plot(
    panels_benz[1], panels_pop[1], panels_pf[1],
    panels_benz[2], panels_pop[2], panels_pf[2],
    panels_benz[3], panels_pop[3], panels_pf[3],
    layout=(3,3), size=(1500, 1200), margin=6Plots.mm)
savefig(fig5_pareto, "h5_pareto_ABC_tijdseries.png")
println("  Saved: h5_pareto_ABC_tijdseries.png")

# ============================================================
#  5b. Glucose-maltose pulse optimisation
# ============================================================
println("\n=== 5b (M3): Glucose-maltose pulse optimisation ===")

function run_gm_puls_m3(N0, t_inf, moi, t_puls, malt_puls, k_tox=0.35)
    p = make_p3b(N0, t_inf, moi, k_tox)
    u0 = zeros(23)
    u0[Model3.Sind] = [4.44, 0.0, 5.42, 0.0]
    u0[Model3.Eind] = [0.95, 0.01, 0.01, 0.01]
    u0[Model3.Nind] = N0
    tspan = (0.0, p.duration)

    maltoseCondition(u,t,integrator)   = t == t_puls
    maltoseAffect!(integrator)         = integrator.u[Model3.Sind[2]] += malt_puls
    maltoseCallBack = DiscreteCallback(maltoseCondition, maltoseAffect!)
    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model3.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)
    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model3.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)
    domainCondition(u,t,integrator)    = any(x->x<0.0, u)
    domainAffect!(integrator)          = Model3.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(Model3.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()), verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([t_puls; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, maltoseCallBack, infectionCallBack, fbaCallBack))
end

dt_waarden   = [0.25, 0.5, 1.0, 1.5, 2.0]
malt_waarden = [0.5, 1.0, 2.337, 5.0]

gm_scores    = zeros(length(dt_waarden), length(malt_waarden))
gm_solutions = Matrix{Any}(nothing, length(dt_waarden), length(malt_waarden))

for (i, dt) in enumerate(dt_waarden)
    for (j, malt) in enumerate(malt_waarden)
        t_puls = opt_t_inf - dt
        t_puls < 0.1 && continue
        sol_gm = run_gm_puls_m3(opt_biomassa, opt_t_inf, opt_moi, t_puls, malt)
        if SciMLBase.successful_retcode(sol_gm)
            B_ts     = [sol_gm.u[k][Model3.Benzind] for k in eachindex(sol_gm.u)]
            B_max    = maximum(B_ts)
            Pf_f     = sol_gm.u[end][Model3.Pfind]
            infectiedosis = opt_moi * opt_biomassa
            idx_95   = findfirst(B_ts .>= 0.95 * B_max)
            t_max_gm = isnothing(idx_95) ? duration : sol_gm.t[idx_95]
            benz_norm_gm = 1.0 - exp(-B_max / BENZ_REF_3b)
            P_norm_gm    = Pf_f > 0.0 ? 1.0 - exp(-log(max(infectiedosis, 2.0)) / log(max(Pf_f, 2.0))) : 0.0
            t_norm_gm    = 1.0 - exp(-max(0.0, t_max_gm - opt_t_inf) / duration)
            gm_scores[i,j]    = 10.0 * benz_norm_gm - 0.1 * P_norm_gm - 0.1 * t_norm_gm
            gm_solutions[i,j] = sol_gm
        end
    end
end

p5b_score = heatmap(string.(malt_waarden), string.(dt_waarden), gm_scores,
    xlabel="Maltose pulse [mmol L⁻¹]",
    ylabel="Δt before infection [h]",
    color=:viridis,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm, right_margin=10Plots.mm)

score_pairs = [(gm_scores[i,j], i, j)
               for i in axes(gm_scores,1), j in axes(gm_scores,2)
               if gm_solutions[i,j] !== nothing]
sort!(score_pairs, by=x->x[1], rev=true)
top4 = score_pairs[1:min(4, length(score_pairs))]

sol_ref_malt = Model3.run(make_p3b(opt_biomassa, opt_t_inf, opt_moi))
B_ref_malt   = [sol_ref_malt.u[i][Model3.Benzind] for i in eachindex(sol_ref_malt.u)]

kleuren_top = [:green :steelblue :darkorange :purple]

p5b_ts = plot(sol_ref_malt.t, B_ref_malt,
    label="Pure maltose (ref)", color=:black, lw=2, linestyle=:dash,
    xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)

for (rank, (sc, i, j)) in enumerate(top4)
    sol_top = gm_solutions[i,j]
    B_top   = [sol_top.u[k][Model3.Benzind] for k in eachindex(sol_top.u)]
    lbl     = "Δt=$(dt_waarden[i])h, malt=$(malt_waarden[j]) (sc=$(round(sc,digits=3)))"
    plot!(p5b_ts, sol_top.t, B_top, label=lbl, color=kleuren_top[rank], lw=2)
end

fig5b_m3 = plot(p5b_score, p5b_ts, layout=(1,2), size=(1200,450), margin=6Plots.mm)
savefig(fig5b_m3, "h5b_m3_gm_puls.png")
println("  Saved: h5b_m3_gm_puls.png")

# ============================================================
#  5c. Toxicity threshold k_tox
# ============================================================
println("\n=== 5c (M3): Toxicity threshold k_tox ===")

k_tox_waarden = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]

benz_ktox = Float64[]; lys_eindktox = Float64[]; t_max_ktox = Float64[]

for kt in k_tox_waarden
    sol_kt = Model3.run(make_p3b(opt_biomassa, opt_t_inf, opt_moi, kt))
    if SciMLBase.successful_retcode(sol_kt)
        B_ts_kt = [sol_kt.u[i][Model3.Benzind] for i in eachindex(sol_kt.u)]
        push!(benz_ktox,    maximum(B_ts_kt))
        push!(lys_eindktox, sol_kt.u[end][Model3.lind])
        push!(t_max_ktox,   sol_kt.t[argmax(B_ts_kt)])
    end
end

pc = plot(k_tox_waarden, benz_ktox, marker=:circle, lw=2, color=:green,
    xlabel="k_tox [h⁻¹]", ylabel="Max Benzonase [mmol L⁻¹]", legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
pd = plot(k_tox_waarden, max.(lys_eindktox, 1.0), marker=:square, lw=2, color=:darkgreen,
    yscale=:log10, xlabel="k_tox [h⁻¹]", ylabel="Lysogenic cells L⁻¹", legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
pe = plot(k_tox_waarden, t_max_ktox, marker=:diamond, lw=2, color=:purple,
    xlabel="k_tox [h⁻¹]", ylabel="Time of max Benzonase [h]", legend=false,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
fig5c_m3 = plot(pc, pd, pe, layout=(1,3), size=(1300,430), margin=6Plots.mm)
savefig(fig5c_m3, "h5c_m3_toxiciteit.png")
println("  Saved: h5c_m3_toxiciteit.png")

# ============================================================
#  5d. Fed-batch with continuous feed rate
# ============================================================
println("\n=== 5d (M3): Fed-batch (continuous feed rate) ===")

C_per_mmol_glc  = 6.0  / 1000.0
C_per_mmol_malt = 12.0 / 1000.0
glc_begin       = 4.44
C_begin         = glc_begin * C_per_mmol_glc

function run_fedbatch_continu_m3(N0, t_inf, moi, feed_rate_mmol_per_h)
    dt_feed  = 1.0/60.0
    amt_feed = feed_rate_mmol_per_h * dt_feed
    p = make_p3b(N0, t_inf, moi)
    u0 = zeros(23)
    u0[Model3.Sind] = [4.44, 0.0, 5.42, 0.0]
    u0[Model3.Eind] = [0.95, 0.01, 0.01, 0.01]
    u0[Model3.Nind] = N0
    tspan = (0.0, p.duration)
    feed_tps = collect(0:dt_feed:p.duration)

    feedCondition(u,t,integrator)      = t in feed_tps
    feedAffect!(integrator)            = integrator.u[Model3.Sind[2]] += amt_feed
    feedCallBack = DiscreteCallback(feedCondition, feedAffect!)
    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model3.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)
    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model3.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)
    domainCondition(u,t,integrator)    = any(x->x<0.0, u)
    domainAffect!(integrator)          = Model3.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    all_tstops = sort(unique(vcat(feed_tps, p.infection_time, fbaUpdateTimepoints)))
    problem = DDEProblem(Model3.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    sol = solve(problem, MethodOfSteps(Tsit5()), verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=all_tstops,
        callback=CallbackSet(domainCallBack, feedCallBack, infectionCallBack, fbaCallBack))

    C_feed   = feed_rate_mmol_per_h * p.duration * C_per_mmol_malt
    C_totaal = C_begin + C_feed
    return sol, C_totaal
end

dosering_waarden = [0.1, 0.5, 1.0, 2.0]
benz_fb = Float64[]; benz_per_cmol = Float64[]
C_batch = 4.44 * C_per_mmol_glc + 2.337 * C_per_mmol_malt
benz_batch_per_cmol = maximum(B_opt) / C_batch

for feed in dosering_waarden
    sol_fb, C_tot = run_fedbatch_continu_m3(opt_biomassa, opt_t_inf, opt_moi, feed)
    B_max = SciMLBase.successful_retcode(sol_fb) ?
        maximum([sol_fb.u[i][Model3.Benzind] for i in eachindex(sol_fb.u)]) : 0.0
    push!(benz_fb, B_max)
    push!(benz_per_cmol, B_max / C_tot)
end

scenario_labels_d = vcat(["Batch\n(ref)"], ["$(d) mmol/L/h" for d in dosering_waarden])
benz_abs_all  = vcat([maximum(B_opt)], benz_fb)
benz_cmol_all = vcat([benz_batch_per_cmol * 1000], benz_per_cmol .* 1000)

p5d_abs = plot(scenario_labels_d, benz_abs_all, marker=:circle, lw=2, color=:teal,
    xlabel="Feed rate [mmol L⁻¹ h⁻¹]", ylabel="Max Benzonase [mmol L⁻¹]",
    legend=false,
    tickfontsize=9, guidefontsize=11,
    bottom_margin=10Plots.mm, left_margin=10Plots.mm)
hline!(p5d_abs, [maximum(B_opt)], color=:black, lw=1.5, linestyle=:dash)

p5d_norm = plot(scenario_labels_d, benz_cmol_all, marker=:square, lw=2, color=:darkorange,
    xlabel="Feed rate [mmol L⁻¹ h⁻¹]", ylabel="Benzonase [µmol C-mol⁻¹]",
    legend=false,
    tickfontsize=9, guidefontsize=11,
    bottom_margin=10Plots.mm, left_margin=10Plots.mm)
hline!(p5d_norm, [benz_batch_per_cmol * 1000], color=:black, lw=1.5, linestyle=:dash)

fig5d_m3 = plot(p5d_abs, p5d_norm, layout=(1,2), size=(1000,450), margin=8Plots.mm)
savefig(fig5d_m3, "h5d_m3_fedbatch.png")
println("  Saved: h5d_m3_fedbatch.png")

# ============================================================
#  5e. Final comparison
# ============================================================
println("\n=== 5e (M3): Final scenario comparison ===")

top1_score, top1_i, top1_j = top4[1]
sol_best_gm  = gm_solutions[top1_i, top1_j]
benz_best_gm = maximum([sol_best_gm.u[i][Model3.Benzind] for i in eachindex(sol_best_gm.u)])
benz_best_fb = length(benz_fb) >= 3 ? benz_fb[3] : benz_fb[end]

scenario_namen = ["Batch", "Opt. batch", "Glc+Malt pulse", "Fed-batch 1 mmol/L/h"]
benz_scenarios = [maximum(B_opt), maximum(B_opt), benz_best_gm, benz_best_fb]

fig5e_m3 = bar(scenario_namen, benz_scenarios,
    color=[:steelblue :green :teal :darkorange], alpha=0.85,
    legend=false,
    ylabel="Max Benzonase [mmol L⁻¹]",
    xlabel="Scenario",
    xrotation=15,
    tickfontsize=11, guidefontsize=13,
    bottom_margin=12Plots.mm, left_margin=10Plots.mm,
    size=(600,430))
savefig(fig5e_m3, "h5e_m3_eindvergelijking.png")

println("\n  Scenario          | Max Benz [mmol/L]")
for i in eachindex(scenario_namen)
    @printf("  %-18s | %.5f\n", scenario_namen[i], benz_scenarios[i])
end
println("  Saved: h5e_m3_eindvergelijking.png")

println("\n=== Chapter 5 Model 3 complete ===")