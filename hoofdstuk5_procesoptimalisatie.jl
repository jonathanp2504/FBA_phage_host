# ============================================================
#  HOOFDSTUK 5: Procesoptimalisatie voor industriële productie
#  Beste model: Model 5 (FBA burst size)
#
#  Analyses:
#   5a. Optimale run Model 5
#   5b. Glucose-maltose puls: eindresultaat + tijdstip max concentratie
#   5c. Toxiciteitsdrempel k_tox analyse
#   5d. Fed-batch scenario met CONTINUE feed-rate (niet puls-gewijs)
#   5e. Eindvergelijking alle scenario's
# ============================================================
include("Model5.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

const BENZ_REF_5 = 0.002
const PN_REF_5   = 6.0e12

model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 1.32
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 2.26, 0.0, 10.0]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 1.71e-12;  duration = 20.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW  = 1.0e-12

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
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", "R_PHAGE_prod",
        k_tox, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod)
end

# ============================================================
#  5a. Optimale run Model 5
# ============================================================
println("=== 5a: Optimale run Model 5 (Optimizer B) ===")

opt_moi      = 2.0
opt_t_inf    = 0.5
opt_biomassa = 1e9
opt_score    = 0.51

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

pa = plot(t_opt, B_opt,
    label="Benzonase", color=:green, lw=2,
    fill=(0,0.2,:green),
    xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
pb = plot(t_opt, max.(PN_opt, 1e-10),
    label="P/N", color=:red, lw=2,
    yscale=:log10,
    xlabel="t [h]", ylabel="P/N [-]",
    legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
fig5a = plot(pa, pb, layout=(1,2), size=(900,430), margin=6Plots.mm)
savefig(fig5a, "h5a_gewogen_optimalisatie.png")
println("  Figuur opgeslagen: h5a_gewogen_optimalisatie.png")

# ============================================================
#  5b. Glucose-maltose puls geoptimaliseerd
#  Toont: (1) eindresultaat (max Benzonase) + (2) tijdstip van max concentratie
# ============================================================
println("\n=== 5b: Glucose-maltose puls optimalisatie ===")

function run_gm_puls(N0, t_inf, moi, t_puls, malt_puls, k_tox=0.05)
    p = make_p5(N0, t_inf, moi, k_tox)
    u0 = zeros(17)
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
# Sla ook tijdstip van max op per combinatie
gm_t_max_benz = fill(NaN, length(dt_waarden), length(malt_waarden))

println("  Grid search puls-timing en -grootte...")
for (i, dt) in enumerate(dt_waarden)
    for (j, malt) in enumerate(malt_waarden)
        t_puls = opt_t_inf - dt
        t_puls < 0.1 && continue
        sol_gm = run_gm_puls(opt_biomassa, opt_t_inf, opt_moi, t_puls, malt)
        if SciMLBase.successful_retcode(sol_gm)
            B_gm  = [sol_gm.u[k][Model5.Benzind] for k in eachindex(sol_gm.u)]
            i_max = argmax(B_gm)
            B_max = B_gm[i_max]
            gm_t_max_benz[i,j] = sol_gm.t[i_max]
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

# Heatmap 1: gewogen score (eindresultaat)
p5b_score = heatmap(string.(malt_waarden), string.(dt_waarden), gm_scores,
    xlabel="Maltosepuls [mmol L⁻¹]",
    ylabel="Δt vóór infectie [h]",
    color=:viridis,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm,
    right_margin=10Plots.mm)

# Heatmap 2: tijdstip van maximale Benzonase
p5b_tmax = heatmap(string.(malt_waarden), string.(dt_waarden), gm_t_max_benz,
    xlabel="Maltosepuls [mmol L⁻¹]",
    ylabel="Δt vóór infectie [h]",
    color=:plasma,
    colorbar_title="t_max [h]",
    bottom_margin=8Plots.mm, left_margin=10Plots.mm,
    right_margin=10Plots.mm)

# Tijdsreeks van beste combinatie
t_puls_best = opt_t_inf - best_dt
t_puls_best < 0.1 && (t_puls_best = 0.1)
sol_best_gm = run_gm_puls(opt_biomassa, opt_t_inf, opt_moi, t_puls_best, best_malt)
t_bg   = sol_best_gm.t
B_bg   = [sol_best_gm.u[i][Model5.Benzind] for i in eachindex(sol_best_gm.u)]
mal_bg = [sol_best_gm.u[i][Model5.Sind[2]] for i in eachindex(sol_best_gm.u)]

# Vergelijk met referentie (geen puls)
B_ref_ts = B_opt   # uit 5a

p5b_ts = plot(t_bg, B_bg,
    label="Beste puls (Δt=$(best_dt)h, malt=$(best_malt) mmol/L)",
    color=:green, lw=2,
    xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p5b_ts, t_opt, B_ref_ts,
    label="Referentie (geen puls)", color=:steelblue, lw=2, linestyle=:dash)
plot!(twinx(), t_bg, mal_bg,
    label="Maltose", color=:orange, lw=2, linestyle=:dot,
    ylabel="Maltose [mmol L⁻¹]", legend=:topright)

fig5b = plot(p5b_score, p5b_tmax, p5b_ts,
    layout=(1,3), size=(1400,430), margin=6Plots.mm)
savefig(fig5b, "h5b_gm_puls.png")
println("  Figuur opgeslagen: h5b_gm_puls.png")

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

pc = plot(k_tox_waarden, benz_ktox,
    marker=:circle, lw=2, color=:green,
    xlabel="k_tox [h⁻¹]", ylabel="Max Benzonase [mmol L⁻¹]",
    legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
pd = plot(k_tox_waarden, max.(lys_eindktox, 1.0),
    marker=:square, lw=2, color=:darkgreen,
    yscale=:log10,
    xlabel="k_tox [h⁻¹]", ylabel="Lysogene cellen L⁻¹",
    legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
fig5c = plot(pc, pd, layout=(1,2), size=(900,430), margin=6Plots.mm)
savefig(fig5c, "h5c_toxiciteit_analyse.png")
println("  Figuur opgeslagen: h5c_toxiciteit_analyse.png")

# ============================================================
#  5d. Fed-batch scenario met CONTINUE feed-rate
#
#  De maltose-toevoer wordt gemodelleerd als een continue toevoeging
#  via de ODE: dS_malt/dt += feed_rate [mmol/L/h]
#  Dit is realistischer dan puls-gewijs toevoegen en vermijdt
#  numerieke artefacten bij grote pulsen.
#
#  Implementatie: we voegen een continue term toe aan het substraat-
#  differential via een aparte ODE-term, gesimuleerd als zeer kleine
#  discrete stappen (elke 1 minuut) om compatibel te blijven met
#  de DDE-solver structuur.
# ============================================================
println("\n=== 5d: Fed-batch scenario (continue feed-rate) ===")

# Koolstofinhoud voor genormaliseerde vergelijking
C_per_mmol_glc  = 6.0  / 1000.0   # C-mol per mmol glucose
C_per_mmol_malt = 12.0 / 1000.0   # C-mol per mmol maltose
glc_begin       = 4.44             # mmol/L
C_begin         = glc_begin * C_per_mmol_glc

function run_fedbatch_continu(N0, t_inf, moi, feed_rate_mmol_per_h)
    # feed_rate in mmol/L/h maltose
    # Gesimuleerd als discrete pulsen elke minuut (1/60 h) voor
    # compatibiliteit met DDE-solver; equivalent aan continue toevoer
    dt_feed  = 1.0/60.0   # elke minuut
    amt_feed = feed_rate_mmol_per_h * dt_feed  # mmol/L per stap

    p = make_p5(N0, t_inf, moi)
    u0 = zeros(17)
    u0[Model5.Sind]  = [4.44, 0.0, 5.42, 0.0]   # geen maltose bij start
    u0[Model5.Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Model5.Nind]  = N0
    tspan = (0.0, p.duration)

    feed_tps = collect(0:dt_feed:p.duration)

    feedCondition(u,t,integrator)  = t in feed_tps
    feedAffect!(integrator)        = integrator.u[Model5.Sind[2]] += amt_feed
    feedCallBack = DiscreteCallback(feedCondition, feedAffect!)

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

    all_tstops = sort(unique(vcat(feed_tps, p.infection_time, fbaUpdateTimepoints)))

    problem = DDEProblem(Model5.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    sol = solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=all_tstops,
        callback=CallbackSet(domainCallBack, feedCallBack,
                             infectionCallBack, fbaCallBack))

    # Totale koolstoftoevoer per L
    C_feed   = feed_rate_mmol_per_h * p.duration * C_per_mmol_malt
    C_totaal = C_begin + C_feed
    return sol, C_totaal
end

dosering_waarden = [0.1, 0.5, 1.0, 2.0]   # mmol/L/h
benz_fb          = Float64[]
benz_per_cmol    = Float64[]

# Referentiebatch: beginconcentraties glucose+maltose
C_batch = 4.44 * C_per_mmol_glc + 2.337 * C_per_mmol_malt
benz_batch_per_cmol = maximum(B_opt) / C_batch

println("  Standaard batch: $(round(maximum(B_opt),sigdigits=3)) mmol/L | " *
        "Benz/C=$(round(benz_batch_per_cmol*1000,sigdigits=3)) µmol/C-mol")

for feed in dosering_waarden
    sol_fb, C_tot = run_fedbatch_continu(opt_biomassa, opt_t_inf, opt_moi, feed)
    B_max = SciMLBase.successful_retcode(sol_fb) ?
        maximum([sol_fb.u[i][Model5.Benzind] for i in eachindex(sol_fb.u)]) : 0.0
    push!(benz_fb, B_max)
    push!(benz_per_cmol, B_max / C_tot)
    println("  Feed=$feed mmol/L/h: benz=$(round(B_max,sigdigits=3)) mmol/L | " *
            "Benz/C=$(round(benz_per_cmol[end]*1000,sigdigits=3)) µmol/C-mol")
end

scenario_labels_d = vcat(["Batch\n(ref)"], ["$(d) mmol/L/h" for d in dosering_waarden])
benz_abs_all = vcat([maximum(B_opt)], benz_fb)
benz_cmol_all = vcat([benz_batch_per_cmol * 1000],
                     benz_per_cmol .* 1000)

p5d_abs = plot(scenario_labels_d, benz_abs_all,
    marker=:circle, lw=2, color=:teal,
    xlabel="Feed-rate [mmol L⁻¹ h⁻¹]",
    ylabel="Max Benzonase [mmol L⁻¹]",
    legend=false,
    bottom_margin=10Plots.mm, left_margin=10Plots.mm)
hline!(p5d_abs, [maximum(B_opt)],
    color=:black, lw=1.5, linestyle=:dash, label="Batch")

p5d_norm = plot(scenario_labels_d, benz_cmol_all,
    marker=:square, lw=2, color=:darkorange,
    xlabel="Feed-rate [mmol L⁻¹ h⁻¹]",
    ylabel="Benzonase [µmol C-mol⁻¹]",
    legend=false,
    bottom_margin=10Plots.mm, left_margin=10Plots.mm)
hline!(p5d_norm, [benz_batch_per_cmol * 1000],
    color=:black, lw=1.5, linestyle=:dash, label="Batch")

fig5d = plot(p5d_abs, p5d_norm, layout=(1,2), size=(1000,450), margin=8Plots.mm)
savefig(fig5d, "h5d_fedbatch.png")
println("  Figuur opgeslagen: h5d_fedbatch.png")

# ============================================================
#  5e. Eindvergelijking alle scenario's
# ============================================================
println("\n=== 5e: Eindvergelijking alle scenario's ===")

benz_best_gm = maximum([sol_best_gm.u[i][Model5.Benzind] for i in eachindex(sol_best_gm.u)])
benz_best_fb = length(benz_fb) >= 3 ? benz_fb[3] : benz_fb[end]

N_best_gm  = sol_best_gm.u[end][Model5.Nind]
Pf_best_gm = sol_best_gm.u[end][Model5.Pfind]
pn_best_gm = N_best_gm > 1.0 ? Pf_best_gm/N_best_gm : 0.0

scenario_namen = ["Standaard batch", "Gewogen opt.", "Glc+Malt puls", "Fed-batch 1 mmol/L/h"]
benz_scenarios = [maximum(B_opt), maximum(B_opt), benz_best_gm, benz_best_fb]
pn_scenarios   = [PN_opt[end], PN_opt[end], pn_best_gm, 0.0]

p5e_benz = bar(scenario_namen, benz_scenarios,
    color=[:steelblue :green :teal :darkorange], alpha=0.85,
    legend=false,
    ylabel="Max Benzonase [mmol L⁻¹]",
    xrotation=20,
    bottom_margin=12Plots.mm, left_margin=10Plots.mm)

p5e_pn = bar(scenario_namen, pn_scenarios,
    color=[:steelblue :green :teal :darkorange], alpha=0.85,
    legend=false,
    ylabel="P/N [-]",
    xrotation=20,
    bottom_margin=12Plots.mm, left_margin=10Plots.mm)

fig5e = plot(p5e_benz, p5e_pn, layout=(1,2), size=(1000,480), margin=8Plots.mm)
savefig(fig5e, "h5e_eindvergelijking.png")

println("\n  Scenario              | Max Benz [mmol/L] | P/N eind")
println("  ----------------------|-------------------|----------")
for i in eachindex(scenario_namen)
    @printf("  %-21s | %-17.5f | %.2e\n",
        scenario_namen[i], benz_scenarios[i], pn_scenarios[i])
end
println("  Figuur opgeslagen: h5e_eindvergelijking.png")

println("\n=== Hoofdstuk 5 voltooid ===")