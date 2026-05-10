# ============================================================
#  HOOFDSTUK 5: Procesoptimalisatie voor industriële productie
#  Beste model: Model 5 (FBA burst size)
#
#  Optimale waarden Model 5 (uit run_optimalisatie_5.jl):
#   Optimizer A: MOI=2.0  | t_inf=0.5h | biomassa=1e9  | score=0.688572
#   Optimizer B: MOI=2.0  | t_inf=0.5h | biomassa=1e7  | score=0.691343
#   → Beste: Optimizer B (hogere score)
#
#  Analyses:
#   5a. Optimale run Model 5
#   5b. Glucose-maltose puls geoptimaliseerd
#   5c. Toxiciteitsdrempel k_tox analyse
#   5d. Fed-batch scenario
#   5e. Eindvergelijking alle scenario's
# ============================================================
include("Model5.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

# Referentiewaarden uit run_kalibratie_5.jl
const BENZ_REF_5 = 0.002
const PN_REF_5   = 6.0e12

# ============================================================
#  Gedeelde parameters
# ============================================================
model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 3.75, 0.0, 4.0]
E_coli_cellDW  = 1e-12
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 6.0e-12;  duration = 20.0
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)

# FBA caches Model 5
naiveModel5   = Model5.loadFBAmodel(model_path)
lysogenModel5 = Model5.addBenzonase!(Model5.loadFBAmodel(model_path), Model5.benz_stoich)
lyticModel5   = Model5.addPhage!(Model5.loadFBAmodel(model_path), Model5.phage_stoich)
all_ex_ids    = [id for id in keys(naiveModel5.reactions) if startswith(id, "R_EX_")]
naiveFba5     = Model5.buildFbaCache(naiveModel5,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba5   = Model5.buildFbaCache(lysogenModel5, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                      benz_id="R_BENZ_prod")
lyticFba5     = Model5.buildLyticFbaCache(lyticModel5, exchange_ids,
                                           "R_BIOMASS_Ec_iJO1366_core_53p95M", "R_PHAGE_prod")

function make_p5(N0, t_inf, moi, k_tox=0.05, f_prod=0.0015)
    Model5.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba5, lysogenFba5, lyticFba5,
        0.0, zeros(4), 0.0, zeros(4), 0.0, 0.0,
        1e-10, 10.0, 5.0, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", "R_PHAGE_prod",
        k_tox, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod)
end

# ============================================================
#  5a. Optimale run Model 5
#  Optimizer B: MOI=2.0 | t_inf=0.5h | biomassa=1e7 | score=0.691343
# ============================================================
println("=== 5a: Optimale run Model 5 (Optimizer B) ===")

opt_moi      = 2.0
opt_t_inf    = 0.5
opt_biomassa = 1e7
opt_score    = 0.691343

println("  MOI=$opt_moi | t_inf=$opt_t_inf h | biomassa=$opt_biomassa | score=$opt_score")

p_opt   = make_p5(opt_biomassa, opt_t_inf, opt_moi)
sol_opt = Model5.run(p_opt)
t_opt   = sol_opt.t
B_opt   = [sol_opt.u[i][Model5.Benzind] for i in eachindex(sol_opt.u)]
Pf_opt  = [sol_opt.u[i][Model5.Pfind]  for i in eachindex(sol_opt.u)]
N_opt   = [sol_opt.u[i][Model5.Nind]   for i in eachindex(sol_opt.u)]
Lys_opt = [sol_opt.u[i][Model5.lind]   for i in eachindex(sol_opt.u)]
PN_opt  = [N_opt[i]>1.0 ? Pf_opt[i]/N_opt[i] : 0.0 for i in eachindex(t_opt)]

println("  Max Benzonase : $(round(maximum(B_opt), sigdigits=3)) mmol/L")
println("  Eind P/N      : $(round(PN_opt[end], sigdigits=2))")

pa = plot(t_opt, B_opt, label="Benzonase", color=:green, lw=2,
    fill=(0,0.2,:green), xlabel="t [h]", ylabel="mmol/L",
    title="5a: Optimale Benzonase (Model 5)")
pb = plot(t_opt, max.(PN_opt, 1e-10), label="P/N", color=:red, lw=2,
    yscale=:log10, xlabel="t [h]", ylabel="P/N [-]",
    title="5a: P/N verhouding (Model 5)")
fig5a = plot(pa, pb, layout=(1,2), size=(900,380))
savefig(fig5a, "h5a_gewogen_optimalisatie.png")
println("  Figuur opgeslagen: h5a_gewogen_optimalisatie.png")

# ============================================================
#  5b. Glucose-maltose puls geoptimaliseerd
# ============================================================
println("\n=== 5b: Glucose-maltose puls optimalisatie ===")

function run_gm_puls(N0, t_inf, moi, t_puls, malt_puls, k_tox=0.05)
    p = make_p5(N0, t_inf, moi, k_tox)
    u0 = zeros(16)
    u0[Model5.Sind]  = [4.44, 0.0, 5.42, 0.0]
    u0[Model5.Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Model5.Nind]  = N0
    tspan = (0.0, p.duration)

    maltoseCondition(u,t,integrator)   = t == t_puls
    maltoseAffect!(integrator)         = integrator.u[Model5.Sind[2]] += malt_puls
    maltoseCallBack = DiscreteCallback(maltoseCondition, maltoseAffect!)

    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model5.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model5.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u,t,integrator)    = any(x->x<0.0, u)
    domainAffect!(integrator)          = Model5.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(Model5.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([t_puls; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, maltoseCallBack,
                             infectionCallBack, fbaCallBack))
end

dt_waarden   = [0.25, 0.5, 1.0, 1.5, 2.0]
malt_waarden = [0.5, 1.0, 2.337, 5.0]

best_gm_score = -Inf
best_dt       = dt_waarden[1]
best_malt     = malt_waarden[1]
gm_scores     = zeros(length(dt_waarden), length(malt_waarden))

println("  Grid search puls-timing en -grootte...")
for (i, dt) in enumerate(dt_waarden)
    for (j, malt) in enumerate(malt_waarden)
        t_puls = opt_t_inf - dt
        t_puls < 0.1 && continue
        sol_gm = run_gm_puls(opt_biomassa, opt_t_inf, opt_moi, t_puls, malt)
        if SciMLBase.successful_retcode(sol_gm)
            B_max = maximum([sol_gm.u[k][Model5.Benzind] for k in eachindex(sol_gm.u)])
            Pf_f  = sol_gm.u[end][Model5.Pfind]
            N_f   = sol_gm.u[end][Model5.Nind]
            PN_f  = N_f > 1.0 ? Pf_f/N_f : 1e12
            sc    = 0.8*(B_max/BENZ_REF_5) - 0.05*(PN_f/PN_REF_5)
            gm_scores[i,j] = sc
            if sc > best_gm_score
                best_gm_score = sc; best_dt = dt; best_malt = malt
            end
        end
    end
end

println("  Beste puls: Δt=$best_dt h voor infectie | malt_puls=$best_malt mmol/L")
println("  Beste score: $(round(best_gm_score, digits=5))")

fig5b_heatmap = heatmap(string.(malt_waarden), string.(dt_waarden), gm_scores,
    xlabel="Maltosepuls [mmol/L]", ylabel="Δt voor infectie [h]",
    title="5b: Gewogen score (puls-timing vs grootte)",
    color=:viridis, size=(550,420))
savefig(fig5b_heatmap, "h5b_gm_puls_heatmap.png")

t_puls_best = opt_t_inf - best_dt
t_puls_best < 0.1 && (t_puls_best = 0.1)
sol_best_gm = run_gm_puls(opt_biomassa, opt_t_inf, opt_moi, t_puls_best, best_malt)
t_bg   = sol_best_gm.t
B_bg   = [sol_best_gm.u[i][Model5.Benzind]  for i in eachindex(sol_best_gm.u)]
mal_bg = [sol_best_gm.u[i][Model5.Sind[2]]  for i in eachindex(sol_best_gm.u)]

fig5b_ts = plot(t_bg, B_bg, label="Benzonase", color=:green, lw=2,
    xlabel="t [h]", ylabel="mmol/L",
    title="5b: Beste puls scenario tijdsreeks")
plot!(twinx(), t_bg, mal_bg, label="Maltose",
    color=:orange, lw=2, linestyle=:dash, ylabel="Maltose [mmol/L]")
savefig(fig5b_ts, "h5b_gm_puls_tijdsreeks.png")
println("  Figuren opgeslagen: h5b_gm_puls_heatmap.png | h5b_gm_puls_tijdsreeks.png")

# ============================================================
#  5c. Toxiciteitsdrempel k_tox analyse
# ============================================================
println("\n=== 5c: Toxiciteitsdrempel k_tox ===")

k_tox_waarden = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]
benz_ktox     = Float64[]
lys_eindktox  = Float64[]

for kt in k_tox_waarden
    sol_kt = Model5.run(make_p5(opt_biomassa, opt_t_inf, opt_moi, kt))
    if SciMLBase.successful_retcode(sol_kt)
        push!(benz_ktox,    maximum([sol_kt.u[i][Model5.Benzind] for i in eachindex(sol_kt.u)]))
        push!(lys_eindktox, sol_kt.u[end][Model5.lind])
        println("  k_tox=$kt: max_benz=$(round(benz_ktox[end],sigdigits=2)) | lys_eind=$(round(lys_eindktox[end],sigdigits=2))")
    end
end

pc = plot(k_tox_waarden, benz_ktox, marker=:circle, lw=2, color=:green,
    xlabel="k_tox [h⁻¹]", ylabel="Max Benzonase [mmol/L]",
    title="5c: Benzonase vs toxiciteitscoëfficiënt", legend=false)
pd = plot(k_tox_waarden, max.(lys_eindktox, 1.0), marker=:square, lw=2, color=:darkgreen,
    yscale=:log10, xlabel="k_tox [h⁻¹]", ylabel="Lysogene cellen/L",
    title="5c: Lysogene populatie vs k_tox", legend=false)
fig5c = plot(pc, pd, layout=(1,2), size=(900,380))
savefig(fig5c, "h5c_toxiciteit_analyse.png")
println("  Figuur opgeslagen: h5c_toxiciteit_analyse.png")

# ============================================================
#  5d. Fed-batch scenario
# ============================================================
println("\n=== 5d: Fed-batch scenario ===")

function run_fedbatch(N0, t_inf, moi, malt_dosering_per_uur)
    p = make_p5(N0, t_inf, moi)
    u0 = zeros(16)
    u0[Model5.Sind]  = [4.44, 0.0, 5.42, 0.0]
    u0[Model5.Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Model5.Nind]  = N0
    tspan = (0.0, p.duration)

    fedbatch_tps = collect(0:0.5:p.duration)
    fedbatch_amt = malt_dosering_per_uur * 0.5

    fedbatchCondition(u,t,integrator)  = t in fedbatch_tps
    fedbatchAffect!(integrator)        = integrator.u[Model5.Sind[2]] += fedbatch_amt
    fedbatchCallBack = DiscreteCallback(fedbatchCondition, fedbatchAffect!)

    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model5.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model5.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u,t,integrator)    = any(x->x<0.0, u)
    domainAffect!(integrator)          = Model5.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(Model5.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([fedbatch_tps; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, fedbatchCallBack,
                             infectionCallBack, fbaCallBack))
end

dosering_waarden = [0.1, 0.5, 1.0, 2.0]
benz_fb = Float64[]
for dos in dosering_waarden
    sol_fb = run_fedbatch(opt_biomassa, opt_t_inf, opt_moi, dos)
    B_fb   = SciMLBase.successful_retcode(sol_fb) ?
        maximum([sol_fb.u[i][Model5.Benzind] for i in eachindex(sol_fb.u)]) : 0.0
    push!(benz_fb, B_fb)
    println("  Dosering=$dos mmol/L/h: max_benz=$(round(B_fb,sigdigits=2)) mmol/L")
end

pe = bar(string.(dosering_waarden), benz_fb,
    xlabel="Maltose dosering [mmol/L/h]", ylabel="Max Benzonase [mmol/L]",
    title="5d: Fed-batch — effect van maltose dosering",
    color=:teal, legend=false)
hline!(pe, [maximum(B_opt)], label="Standaard (batch)",
    color=:black, lw=2, linestyle=:dash)
savefig(pe, "h5d_fedbatch.png")
println("  Figuur opgeslagen: h5d_fedbatch.png")

# ============================================================
#  5e. Eindvergelijking alle scenario's
# ============================================================
println("\n=== 5e: Eindvergelijking alle scenario's ===")

benz_best_gm  = maximum([sol_best_gm.u[i][Model5.Benzind] for i in eachindex(sol_best_gm.u)])
benz_best_fb  = length(benz_fb) >= 3 ? benz_fb[3] : benz_fb[end]
N_best_gm     = sol_best_gm.u[end][Model5.Nind]
Pf_best_gm    = sol_best_gm.u[end][Model5.Pfind]
pn_best_gm    = N_best_gm > 1.0 ? Pf_best_gm/N_best_gm : 0.0

scenario_namen = ["Standaard\n(batch)" "Gewogen\noptimalisatie" "Glc+Malt puls" "Fed-batch\n(1 mmol/L/h)"]
benz_scenarios = [maximum(B_opt), maximum(B_opt), benz_best_gm, benz_best_fb]
pn_scenarios   = [PN_opt[end], PN_opt[end], pn_best_gm, 0.0]

pn_vals = filter(x -> !isinf(x) && x > 0, pn_scenarios)
pn_max  = isempty(pn_vals) ? 1.0 : maximum(pn_vals)

scenario_labels = ["Standaard", "Gewogen opt.", "Glc+Malt puls", "Fed-batch (1 mmol/h)"]

fig5e = bar(scenario_labels, benz_scenarios,
    color=[:steelblue :green :teal :darkorange],
    title="5e: Eindvergelijking scenario's (Model 5)",
    ylabel="Max Benzonase [mmol/L]",
    xlabel="Scenario",
    legend=false, size=(750,430))
savefig(fig5e, "h5e_eindvergelijking.png")

println("\n  Scenario              | Max Benz [mmol/L] | P/N eind")
println("  ----------------------|-------------------|----------")
for i in eachindex(scenario_namen)
    @printf("  %-21s | %-17.5f | %.2e\n",
        replace(scenario_namen[i], "\n"=>" "),
        benz_scenarios[i], pn_scenarios[i])
end

println("\n=== Hoofdstuk 5 voltooid ===")