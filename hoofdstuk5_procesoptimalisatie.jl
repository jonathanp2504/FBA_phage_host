# ============================================================
#  HOOFDSTUK 5: Procesoptimalisatie voor industriële productie
#  Beste model uit Hoofdstuk 4 + glucose-maltose puls
#
#  Analyses:
#   5a. Volledige gewogen optimalisatie (best model)
#   5b. Glucose-maltose puls geoptimaliseerd (t_puls + malt_puls)
#   5c. Toxiciteitsdrempel k_tox analyse
#   5d. Fed-batch scenario: continue maltose dosering
#   5e. Eindvergelijking: alle scenario's naast elkaar
# ============================================================

using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels, OrdinaryDiffEqCore

if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end

# Laad het beste model (Model 3 als uitgangspunt, pas aan naar best uit H4)
include("./Bin_model_3/parameters.jl")
include("./Bin_model_3/FBA.jl")
include("./Bin_model_3/dFBA.jl")
include("optimization_weighted_B.jl")

model_path    = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn     = 2.0;  beta_deg = 0.5
K_s           = [0.0278, 0.0146, 0.0543, 0.0833]
tau           = 1.0;  b = 170.0
p_pref        = [0.8925, 0.08925, 0.008925, 0.008925]
V_max         = [0.0, 3.75, 0.0, 4.0]
E_coli_cellDW = 1e-12
exchange_ids  = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values     = [180.16, 342.3, 92.09, 60.05]
h_release     = 6.0e-12;  duration = 20.0
mu_max_vector = [1.33, 1.26, 1.10, 0.29]
e_max_vector  = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)

naiveModel   = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
all_ex_ids   = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
naiveFba     = buildFbaCache(naiveModel,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba   = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                              benz_id="R_BENZ_prod")

function make_p(N0, t_inf, moi, k_tox=0.05, f_prod=0.0015)
    Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba, lysogenFba, 0.0, zeros(4), 0.0, zeros(4), 0.0,
        1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
        "R_BENZ_prod", k_tox, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)
end

# ============================================================
#  5a. Volledige gewogen optimalisatie
# ============================================================
println("=== 5a: Gewogen optimalisatie (80% Benzonase, 20% P/N) ===")
p_base = make_p(1e9, 2.0, 2.0)
opt    = run_optimization_B(p_base)

p_opt   = make_p(opt.best_biomass, opt.best_t_inf, FIXED_MOI_B)
sol_opt = run(p_opt)
t_opt   = sol_opt.t
B_opt   = [sol_opt.u[i][Benzind] for i in eachindex(sol_opt.u)]
Pf_opt  = [sol_opt.u[i][Pfind]  for i in eachindex(sol_opt.u)]
N_opt   = [sol_opt.u[i][Nind]   for i in eachindex(sol_opt.u)]
Lys_opt = [sol_opt.u[i][lind]   for i in eachindex(sol_opt.u)]
PN_opt  = [N_opt[i]>1.0 ? Pf_opt[i]/N_opt[i] : 0.0 for i in eachindex(t_opt)]

println("  Max Benzonase : $(round(maximum(B_opt), sigdigits=3)) mmol/L")
println("  Eind P/N      : $(round(PN_opt[end], sigdigits=2))")

pa = plot(t_opt, B_opt, label="Benzonase", color=:green, lw=2,
    fill=(0,0.2,:green), xlabel="t [h]", ylabel="mmol/L",
    title="5a: Optimale Benzonase (gewogen obj.)")
pb = plot(t_opt, PN_opt, label="P/N", color=:red, lw=2,
    yscale=:log10, xlabel="t [h]", ylabel="P/N [-]",
    title="5a: P/N verhouding (gewogen obj.)")
fig5a = plot(pa, pb, layout=(1,2), size=(900,380))
savefig(fig5a, "h5a_gewogen_optimalisatie.png")
println("  Figuur opgeslagen: h5a_gewogen_optimalisatie.png")

# ============================================================
#  5b. Glucose-maltose puls geoptimaliseerd
#  Extra parameters: t_puls (relatief tot t_inf) en malt_puls
# ============================================================
println("\n=== 5b: Glucose-maltose puls optimalisatie ===")

function run_gm_puls(N0, t_inf, moi, t_puls, malt_puls, k_tox=0.05)
    p = make_p(N0, t_inf, moi, k_tox)
    u0 = zeros(23)
    u0[Sind]  = [4.44, 0.0, 5.42, 0.0]   # glucose AAN, maltose UIT
    u0[Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Nind]  = N0
    tspan     = (0.0, p.duration)

    maltoseCondition(u, t, integrator)   = t == t_puls
    maltoseAffect!(integrator)           = integrator.u[Sind[2]] += malt_puls
    maltoseCallBack = DiscreteCallback(maltoseCondition, maltoseAffect!)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator)  = t in fbaUpdateTimepoints
    fbaAffect!(integrator)              = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u,t,integrator)     = any(x->x<0.0, u)
    domainAffect!(integrator)           = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([t_puls; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, maltoseCallBack,
                             infectionCallBack, fbaCallBack))
end

# Sweep: tijdstip puls (0.5 tot 2.0 uur voor infectie) en grootte puls
t_inf_gm      = opt.best_t_inf
N0_gm         = opt.best_biomass
dt_waarden    = [0.25, 0.5, 1.0, 1.5, 2.0]
malt_waarden  = [0.5, 1.0, 2.337, 5.0]

best_gm_score = -Inf
best_dt       = dt_waarden[1]
best_malt     = malt_waarden[1]
gm_scores     = zeros(length(dt_waarden), length(malt_waarden))

println("  Grid search puls-timing en -grootte...")
for (i, dt) in enumerate(dt_waarden)
    for (j, malt) in enumerate(malt_waarden)
        t_puls = t_inf_gm - dt
        t_puls < 0.1 && continue
        sol_gm = run_gm_puls(N0_gm, t_inf_gm, 2.0, t_puls, malt)
        if SciMLBase.successful_retcode(sol_gm)
            B_max = maximum([sol_gm.u[k][Benzind] for k in eachindex(sol_gm.u)])
            Pf_f  = sol_gm.u[end][Pfind]; N_f = sol_gm.u[end][Nind]
            PN_f  = N_f>1.0 ? Pf_f/N_f : 1e12
            sc    = 0.8*(B_max/BENZ_REF_B) - 0.2*(PN_f/PN_REF_B)
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

# Tijdsreeks beste puls scenario
sol_best_gm = run_gm_puls(N0_gm, t_inf_gm, 2.0, t_inf_gm-best_dt, best_malt)
t_bg  = sol_best_gm.t
B_bg  = [sol_best_gm.u[i][Benzind] for i in eachindex(sol_best_gm.u)]
mal_bg = [sol_best_gm.u[i][Sind[2]] for i in eachindex(sol_best_gm.u)]

fig5b_ts = plot(t_bg, B_bg, label="Benzonase", color=:green, lw=2,
    xlabel="t [h]", ylabel="mmol/L",
    title="5b: Beste puls scenario tijdsreeks")
plot!(twinx(), t_bg, mal_bg, label="Maltose",
    color=:orange, lw=2, linestyle=:dash, ylabel="Maltose [mmol/L]")
savefig(fig5b_ts, "h5b_gm_puls_tijdsreeks.png")
println("  Figuren opgeslagen: h5b_gm_puls_heatmap.png | h5b_gm_puls_tijdsreeks.png")

# ============================================================
#  5c. Toxiciteitsdrempel analyse
# ============================================================
println("\n=== 5c: Toxiciteitsdrempel k_tox ===")

k_tox_waarden = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]
benz_ktox     = Float64[]
lys_eindktox  = Float64[]

for kt in k_tox_waarden
    sol_kt = run(make_p(opt.best_biomass, opt.best_t_inf, 2.0, kt))
    if SciMLBase.successful_retcode(sol_kt)
        push!(benz_ktox,    maximum([sol_kt.u[i][Benzind] for i in eachindex(sol_kt.u)]))
        push!(lys_eindktox, sol_kt.u[end][lind])
        println("  k_tox=$kt: max_benz=$(round(benz_ktox[end],sigdigits=2)) | lys_eind=$(round(lys_eindktox[end],sigdigits=2))")
    end
end

pc = plot(k_tox_waarden, benz_ktox,
    marker=:circle, lw=2, color=:green,
    xlabel="k_tox [h⁻¹]", ylabel="Max Benzonase [mmol/L]",
    title="5c: Benzonase vs toxiciteitscoëfficiënt", legend=false)
pd = plot(k_tox_waarden, lys_eindktox,
    marker=:square, lw=2, color=:darkgreen,
    yscale=:log10, xlabel="k_tox [h⁻¹]", ylabel="Lysogene cellen/L",
    title="5c: Lysogene populatie vs k_tox", legend=false)
fig5c = plot(pc, pd, layout=(1,2), size=(900,380))
savefig(fig5c, "h5c_toxiciteit_analyse.png")
println("  Figuur opgeslagen: h5c_toxiciteit_analyse.png")

# ============================================================
#  5d. Fed-batch: continue lage maltosedosering
# ============================================================
println("\n=== 5d: Fed-batch scenario ===")

function run_fedbatch(N0, t_inf, moi, malt_dosering_per_uur)
    p = make_p(N0, t_inf, moi)
    u0 = zeros(23)
    u0[Sind]  = [4.44, 0.0, 5.42, 0.0]
    u0[Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Nind]  = N0
    tspan     = (0.0, p.duration)

    # Voeg elke 0.5 uur een kleine maltosedosis toe vanaf t=0
    fedbatch_tps = collect(0:0.5:p.duration)
    fedbatch_amt = malt_dosering_per_uur * 0.5

    fedbatchCondition(u,t,integrator)  = t in fedbatch_tps
    fedbatchAffect!(integrator)        = integrator.u[Sind[2]] += fedbatch_amt
    fedbatchCallBack = DiscreteCallback(fedbatchCondition, fedbatchAffect!)

    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u,t,integrator)    = any(x->x<0.0, u)
    domainAffect!(integrator)          = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([fedbatch_tps; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, fedbatchCallBack,
                             infectionCallBack, fbaCallBack))
end

dosering_waarden = [0.1, 0.5, 1.0, 2.0]
benz_fb = Float64[]
for dos in dosering_waarden
    sol_fb = run_fedbatch(opt.best_biomass, opt.best_t_inf, 2.0, dos)
    B_fb   = SciMLBase.successful_retcode(sol_fb) ?
        maximum([sol_fb.u[i][Benzind] for i in eachindex(sol_fb.u)]) : 0.0
    push!(benz_fb, B_fb)
    println("  Dosering=$dos mmol/L/h: max_benz=$(round(B_fb,sigdigits=2)) mmol/L")
end

pe = bar(string.(dosering_waarden), benz_fb,
    xlabel="Maltose dosering [mmol/L/h]",
    ylabel="Max Benzonase [mmol/L]",
    title="5d: Fed-batch — effect van maltose dosering",
    color=:teal, legend=false)
hline!(pe, [maximum(B_opt)], label="Standaard (maltose als batch)",
    color=:black, lw=2, linestyle=:dash)
savefig(pe, "h5d_fedbatch.png")
println("  Figuur opgeslagen: h5d_fedbatch.png")

# ============================================================
#  5e. Eindvergelijking alle scenario's
# ============================================================
println("\n=== 5e: Eindvergelijking alle scenario's ===")

scenario_namen  = ["Standaard\n(batch maltose)" "Gewogen\noptimalisatie"
                   "Glc+Malt puls" "Fed-batch\n(1 mmol/L/h)"]
benz_scenarios  = [
    maximum(B3),                          # standaard Model 3
    maximum(B_opt),                       # gewogen optimalisatie
    maximum([sol_best_gm.u[i][Benzind] for i in eachindex(sol_best_gm.u)]),
    length(benz_fb) >= 3 ? benz_fb[3] : benz_fb[end]
]
pn_scenarios = [
    PN_opt[end],
    PN_opt[end],
    N_opt[end] > 1.0 ? sol_best_gm.u[end][Pfind]/sol_best_gm.u[end][Nind] : 0.0,
    0.0
]

fig5e = groupedbar(
    scenario_namen,
    hcat(benz_scenarios, pn_scenarios ./ maximum(filter(!isinf, pn_scenarios)) .* maximum(benz_scenarios))',
    label=["Max Benzonase [mmol/L]" "P/N (genorm. op schaal Benz)"],
    color=[:green :red], alpha=0.8,
    title="5e: Eindvergelijking scenario's",
    ylabel="Waarde", legend=:topright, size=(750,430))
savefig(fig5e, "h5e_eindvergelijking.png")

println("\n  Scenario              | Max Benz [mmol/L] | P/N eind")
println("  ----------------------|-------------------|----------")
for i in eachindex(scenario_namen)
    @printf("  %-21s | %-17.5f | %.2e\n",
        replace(scenario_namen[i], "\n"=>" "),
        benz_scenarios[i], pn_scenarios[i])
end

println("\n=== Hoofdstuk 5 voltooid ===")
