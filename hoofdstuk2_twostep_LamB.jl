# ============================================================
#  HOOFDSTUK 2: Two-step adsorptie en de rol van LamB
#
#  Analyses:
#   2a. Populatiedynamica Model 2 bij 4 MOI-waarden
#   2b. NIEUW: Model 1 met alfa van two-step — toont waarom one-step
#       faalt bij hoge concentratie (injectie rate-limiting)
#   2c. Glucose-maltose puls: eindresultaat + tijdstip max concentratie
#   2d. Lage beginbiomassa simulatie
#   2e. Fagen voor maltosepuls (delta_t sweep): eindresultaat + tijdstip
#   2f. Benzonase productie bij 4 MOI-waarden
# ============================================================
include("Model1.jl")
include("Model2.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

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
h_release      = 1.71e-12;  duration = 20.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW  = 1.0e-12

baseModel2  = Model2.loadFBAmodel(model_path)
all_ex_ids  = [id for id in keys(baseModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2   = Model2.buildFbaCache(baseModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2 = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                    "R_BIOMASS_Ec_iJO1366_core_53p95M")
fbaCache1   = Model1.buildFbaCache(Model1.loadFBAmodel(model_path), exchange_ids,
                                    "R_BIOMASS_Ec_iJO1366_core_53p95M")

function make_p2(N0, t_inf, moi)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0,
        mu_max_vector, e_max_vector, 0.0015, 0.05)
end

function make_p1(N0, t_inf, moi; alfa_=1.3e-9)
    Model1.Parameters(duration, N0, alpha_syn, beta_deg, K_s,
        [7.24, 2.26, 0.0, 10.0], p_pref, tau, b, alfa_, 1.0e-12,
        MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        fbaCache1, 0.0, zeros(4), t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0, mu_max_vector, e_max_vector, 0.0, 0.0)
end

# ============================================================
#  2a. Populatiedynamica Model 2 bij 4 MOI-waarden
# ============================================================
println("=== 2a: Model 2 populatiedynamica bij verschillende MOI ===")

moi_compare = [0.001, 0.1, 1.0, 5.0]
panels = []

for moi in moi_compare
    sol2 = Model2.run(make_p2(1e9, 2.0, moi))
    t2   = sol2.t
    N2   = [sol2.u[i][Model2.Nind]  for i in eachindex(sol2.u)]
    Lyt2 = [sol2.u[i][Model2.Lind]  for i in eachindex(sol2.u)]
    Lys2 = [sol2.u[i][Model2.lind]  for i in eachindex(sol2.u)]
    Pf2  = [max(sol2.u[i][Model2.Pfind], 1.0) for i in eachindex(sol2.u)]
    Pa2  = [sol2.u[i][Model2.Paind] for i in eachindex(sol2.u)]

    p_panel = plot(t2, [N2 Lyt2 Lys2 Pf2],
        label=["Naïef N" "Lytisch L" "Lysogeen l" "Vrije fagen Pf"],
        color=[:blue :red :green :darkred],
        lw=2, yscale=:log10, ylims=(1,:auto),
        xlabel="t [h]", ylabel="Cellen of fagen L⁻¹",
        legend = moi==moi_compare[1] ? :topright : false,
        bottom_margin=5Plots.mm, left_margin=8Plots.mm)
    plot!(p_panel, t2, max.(Pa2, 1.0), label="Gehecht Pa",
        color=:purple, lw=1.5, linestyle=:dot,
        legend = moi==moi_compare[1] ? :topright : false)
    push!(panels, p_panel)
end

fig2a = plot(panels..., layout=(2,2), size=(1000,720), margin=5Plots.mm)
savefig(fig2a, "h2a_model2_moi_vergelijking.png")
println("  Figuur opgeslagen: h2a_model2_moi_vergelijking.png")

# ============================================================
#  2b. NIEUW: One-step model met alfa van two-step
#
#  Bij LAGE MOI is diffusie naar het celoppervlak de snelheidsbeperkende
#  stap: een lage adsorptieconstante weerspiegelt de kans op een productieve
#  botsing en het one-step model is een redelijke benadering.
#
#  Bij HOGE concentraties (hoge k_att zoals in het two-step model) is
#  de fysische aanhechtingsstap snel: de injectie zelf wordt dan de
#  snelheidsbeperkende stap. Het one-step model kan dit niet vatten —
#  het instantaan infecteren van de hele populatie veroorzaakt numerieke
#  instabiliteit of biologisch onrealistische dynamica.
#
#  Test: draai Model 1 met alfa_ads = k_att van Model 2 (7.92e-8 L/cel/h)
#  en vergelijk met Model 2 bij dezelfde parameter.
# ============================================================
println("\n=== 2b: One-step met alfa van two-step vs Model 2 ===")

alfa_twostep = 7.92e-8   # k_att van Model 2
moi_test     = 2.0       # MOI waarbij het verschil het grootst is

println("  Model 1 (one-step) met alfa_ads=$(alfa_twostep)...")
sol_m1_hoog = Model1.run(make_p1(1e9, 2.0, moi_test; alfa_=alfa_twostep))
println("  Model 2 (two-step) met k_att=$(alfa_twostep)...")
sol_m2_ref  = Model2.run(make_p2(1e9, 2.0, moi_test))
println("  Model 1 (one-step) met originele alfa_ads=1.3e-9...")
sol_m1_laag = Model1.run(make_p1(1e9, 2.0, moi_test; alfa_=1.3e-9))

# Hulpfunctie: extraheer populatievariabelen uit Model 1-oplossing
function extract_m1(sol)
    t  = sol.t
    S  = [sol.u[i][Model1.Sind_S] for i in eachindex(sol.u)]
    I  = [sol.u[i][Model1.Sind_I] for i in eachindex(sol.u)]
    L  = [sol.u[i][Model1.Sind_L] for i in eachindex(sol.u)]
    P  = max.([sol.u[i][Model1.Pfind] for i in eachindex(sol.u)], 1.0)
    X  = S .+ I .+ L
    return t, S, I, L, P, X
end

function extract_m2(sol)
    t   = sol.t
    N   = [sol.u[i][Model2.Nind]  for i in eachindex(sol.u)]
    Lyt = [sol.u[i][Model2.Lind]  for i in eachindex(sol.u)]
    Lys = [sol.u[i][Model2.lind]  for i in eachindex(sol.u)]
    P   = max.([sol.u[i][Model2.Pfind] for i in eachindex(sol.u)], 1.0)
    X   = N .+ Lyt .+ Lys
    return t, N, Lyt, Lys, P, X
end

t1h, S1h, I1h, L1h, P1h, X1h = extract_m1(sol_m1_hoog)
t2,  N2,  Lyt2, Lys2, P2, X2 = extract_m2(sol_m2_ref)
t1l, S1l, I1l, L1l, P1l, X1l = extract_m1(sol_m1_laag)

# Panel 1: totale biomassa
p2b_X = plot(xlabel="t [h]", ylabel="Totale biomassa (cellen L⁻¹)",
    yscale=:log10, ylims=(1e4,:auto), legend=:bottomleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p2b_X, t1l, max.(X1l, 1.0),
    label="Model 1 (α=1.3×10⁻⁹)", color=:steelblue,  lw=2, linestyle=:dash)
plot!(p2b_X, t1h, max.(X1h, 1.0),
    label="Model 1 (α=k_att two-step)", color=:red,    lw=2, linestyle=:dot)
plot!(p2b_X, t2,  max.(X2,  1.0),
    label="Model 2 (two-step)",  color=:darkgreen, lw=2)

# Panel 2: vrije fagen
p2b_P = plot(xlabel="t [h]", ylabel="Vrije fagen L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p2b_P, t1l, P1l,
    label="Model 1 (α=1.3×10⁻⁹)", color=:steelblue,  lw=2, linestyle=:dash)
plot!(p2b_P, t1h, P1h,
    label="Model 1 (α=k_att two-step)", color=:red,    lw=2, linestyle=:dot)
plot!(p2b_P, t2,  P2,
    label="Model 2 (two-step)",  color=:darkgreen, lw=2)

# Panel 3: lysogene fractie
function lysogene_frac_m1(sol)
    frac = Float64[]
    for i in eachindex(sol.u)
        I_ = sol.u[i][Model1.Sind_I]
        L_ = sol.u[i][Model1.Sind_L]
        tot = I_ + L_
        push!(frac, tot > 1.0 ? L_/tot : 0.0)
    end
    return frac
end
function lysogene_frac_m2(sol)
    frac = Float64[]
    for i in eachindex(sol.u)
        Lyt_ = sol.u[i][Model2.Lind]
        Lys_ = sol.u[i][Model2.lind]
        tot  = Lyt_ + Lys_
        push!(frac, tot > 1.0 ? Lys_/tot : 0.0)
    end
    return frac
end

p2b_frac = plot(xlabel="t [h]", ylabel="Lysogene fractie [-]",
    ylims=(0,1.05), legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p2b_frac, t1l, lysogene_frac_m1(sol_m1_laag),
    label="Model 1 (α=1.3×10⁻⁹)", color=:steelblue,  lw=2, linestyle=:dash)
plot!(p2b_frac, t1h, lysogene_frac_m1(sol_m1_hoog),
    label="Model 1 (α=k_att two-step)", color=:red,    lw=2, linestyle=:dot)
plot!(p2b_frac, t2,  lysogene_frac_m2(sol_m2_ref),
    label="Model 2 (two-step)",  color=:darkgreen, lw=2)

fig2b = plot(p2b_X, p2b_P, p2b_frac, layout=(1,3), size=(1300,430), margin=6Plots.mm)
savefig(fig2b, "h2b_onestep_vs_twostep_alfa.png")
println("  Figuur opgeslagen: h2b_onestep_vs_twostep_alfa.png")

# ============================================================
#  2c. Glucose-maltose puls scenario (Model 2)
#  Toont: (1) eindresultaat Benzonase + (2) tijdstip van max concentratie
# ============================================================
println("\n=== 2c: Glucose-maltose puls scenario ===")

function run_glucose_maltose(N0, t_inf, moi, t_puls, malt_puls)
    p = make_p2(N0, t_inf, moi)
    u0 = zeros(23)
    u0[Model2.Sind]  = [4.44, 0.0, 5.42, 0.0]
    u0[Model2.Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Model2.Nind]  = N0
    tspan = (0.0, p.duration)

    maltoseCondition(u,t,integrator)   = t == t_puls
    maltoseAffect!(integrator)         = integrator.u[Model2.Sind[2]] += malt_puls
    maltoseCallBack = DiscreteCallback(maltoseCondition, maltoseAffect!)

    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model2.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model2.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u,t,integrator) = any(x->x<0.0, u)
    domainAffect!(integrator)       = Model2.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(Model2.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([t_puls; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, maltoseCallBack, infectionCallBack, fbaCallBack))
end

t_inf_glc      = 5.0
malt_puls      = 2.337
delta_t_values = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
max_benz_glc   = Float64[]   # eindresultaat Benzonase
t_max_benz_glc = Float64[]   # tijdstip van maximale Benzonase

println("  Maltosepuls sweep (t_inf=$t_inf_glc h, MOI=2.0)")
for dt in delta_t_values
    t_puls = t_inf_glc - dt
    t_puls < 0.1 && continue
    sol_gm = run_glucose_maltose(1e9, t_inf_glc, 2.0, t_puls, malt_puls)
    if SciMLBase.successful_retcode(sol_gm)
        B_ts   = [sol_gm.u[i][Model2.Benzind] for i in eachindex(sol_gm.u)]
        i_max  = argmax(B_ts)
        push!(max_benz_glc,   B_ts[i_max])
        push!(t_max_benz_glc, sol_gm.t[i_max])
        println("  dt=$dt h: max_benz=$(round(B_ts[i_max],sigdigits=3)) mmol/L | ",
                "t_max=$(round(sol_gm.t[i_max],digits=2)) h")
    else
        push!(max_benz_glc,   0.0)
        push!(t_max_benz_glc, NaN)
        println("  dt=$dt h: simulatie gefaald")
    end
end

# Referentie: altijd maltose aanwezig
sol_std   = Model2.run(make_p2(1e9, t_inf_glc, 2.0))
B_std     = [sol_std.u[i][Model2.Benzind] for i in eachindex(sol_std.u)]
Benz_std  = maximum(B_std)
t_std_max = sol_std.t[argmax(B_std)]

dt_labels = string.(delta_t_values[1:length(max_benz_glc)])

# Panel A: eindresultaat (max Benzonase)
p2c_eind = plot(dt_labels, max_benz_glc,
    seriestype=:bar, color=:teal, legend=:topright,
    xlabel="Δt maltosepuls vóór infectie [h]",
    ylabel="Max Benzonase [mmol L⁻¹]",
    label="Glc+Malt puls",
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
hline!(p2c_eind, [Benz_std],
    color=:black, lw=2, linestyle=:dash, label="Altijd maltose (ref)")

# Panel B: tijdstip van max Benzonase
p2c_tijdstip = plot(dt_labels, t_max_benz_glc,
    seriestype=:bar, color=:darkorange, legend=:topright,
    xlabel="Δt maltosepuls vóór infectie [h]",
    ylabel="Tijdstip max Benzonase [h]",
    label="Glc+Malt puls",
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
hline!(p2c_tijdstip, [t_std_max],
    color=:black, lw=2, linestyle=:dash, label="Altijd maltose (ref)")

# Panel C: tijdsreeks beste puls
idx_best = argmax(max_benz_glc)
dt_best  = delta_t_values[idx_best]
sol_best = run_glucose_maltose(1e9, t_inf_glc, 2.0, t_inf_glc-dt_best, malt_puls)
t_b      = sol_best.t
B_b      = [sol_best.u[i][Model2.Benzind] for i in eachindex(sol_best.u)]
N_b      = [sol_best.u[i][Model2.Nind]    for i in eachindex(sol_best.u)]
Lys_b    = [sol_best.u[i][Model2.lind]    for i in eachindex(sol_best.u)]
mal_b    = [sol_best.u[i][Model2.Sind[2]] for i in eachindex(sol_best.u)]

p2c_ts = plot(t_b, B_b,
    label="Benzonase [mmol/L]", color=:green, lw=2,
    xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(twinx(), t_b, mal_b,
    label="Maltose [mmol/L]", color=:orange, lw=2, linestyle=:dash,
    ylabel="Maltose [mmol L⁻¹]", legend=:topright)

fig2c = plot(p2c_eind, p2c_tijdstip, p2c_ts,
    layout=(1,3), size=(1300,430), margin=6Plots.mm)
savefig(fig2c, "h2c_glucose_maltose_puls.png")
println("  Figuur opgeslagen: h2c_glucose_maltose_puls.png")

# ============================================================
#  2d. Lage beginbiomassa simulatie (N0=1e3, infectie bij N≈1e10)
# ============================================================
println("\n=== 2d: Lage beginbiomassa (N0=1e3, infectie bij N≈1e10) ===")

p_groei = make_p2(1e3, 999.0, 0.0)
u0_groei = zeros(23)
u0_groei[Model2.Sind] = [4.44, 2.337, 5.42, 0.0]
u0_groei[Model2.Eind] = [0.95, 0.01, 0.01, 0.01]
u0_groei[Model2.Nind] = 1e3

fbaUpdateTimepoints_g = collect(0:1/60:15.0)
fbaUpdateCondition_g(u,t,integrator) = t in fbaUpdateTimepoints_g
fbaAffect_g!(integrator)             = Model2.fbaUpdate!(integrator.u, p_groei)
fbaCallBack_g = DiscreteCallback(fbaUpdateCondition_g, fbaAffect_g!)
domainCondition_g(u,t,integrator)    = any(x->x<0.0, u)
domainAffect_g!(integrator)          = Model2.enforcePositiveDomain!(integrator.u)
domainCallBack_g = DiscreteCallback(domainCondition_g, domainAffect_g!)

prob_groei = DDEProblem(Model2.simulate_dFBA!, u0_groei, (p_groei,t)->u0_groei,
                        (0.0, 15.0), p_groei)
sol_groei = solve(prob_groei, MethodOfSteps(Tsit5()),
    verbose=false, reltol=1e-4, abstol=1e-6,
    tstops=fbaUpdateTimepoints_g,
    callback=CallbackSet(domainCallBack_g, fbaCallBack_g))

N_groei = [sol_groei.u[i][Model2.Nind] for i in eachindex(sol_groei.u)]
idx_1e10 = findfirst(x -> x >= 1e10, N_groei)
t_infectie_laat = isnothing(idx_1e10) ? sol_groei.t[end] : sol_groei.t[idx_1e10]
println("  Biomassa bereikt 1e10 op t=$(round(t_infectie_laat, digits=2)) h")

p_laat = make_p2(1e3, t_infectie_laat, 2.0)
p_laat_lang = Model2.Parameters(t_infectie_laat + 15.0, p_laat.startingBiomass,
    p_laat.alpha_syn, p_laat.beta_deg, p_laat.K_s, p_laat.V_max, p_laat.p_pref,
    p_laat.tau, p_laat.b, p_laat.E_coli_cellDW, p_laat.MW, p_laat.h_release,
    p_laat.biomass_id, p_laat.ex_ids, p_laat.all_exchanges, p_laat.essentials,
    p_laat.fbaModelNaive, p_laat.fbaModelLysogen,
    p_laat.mu_N, copy(p_laat.q_N), p_laat.mu_l, copy(p_laat.q_l), p_laat.q_benz_l,
    p_laat.k_attach, p_laat.k_dettach, p_laat.k_inject, p_laat.K_mal,
    t_infectie_laat, 2.0 * 1e3,
    p_laat.benz_id, p_laat.k_tox, p_laat.beta_benz, p_laat.q_benz,
    copy(p_laat.mu_max), copy(p_laat.e_max), p_laat.f_prod, p_laat.Y_benz)

u0_laat = zeros(23)
u0_laat[Model2.Sind] = [4.44, 2.337, 5.42, 0.0]
u0_laat[Model2.Eind] = [0.95, 0.01, 0.01, 0.01]
u0_laat[Model2.Nind] = 1e3
sol_laat = Model2.run(p_laat_lang)

if SciMLBase.successful_retcode(sol_laat)
    t_ld  = sol_laat.t
    N_ld  = [sol_laat.u[i][Model2.Nind]    for i in eachindex(sol_laat.u)]
    Lys_ld= [sol_laat.u[i][Model2.lind]    for i in eachindex(sol_laat.u)]
    Lyt_ld= [sol_laat.u[i][Model2.Lind]    for i in eachindex(sol_laat.u)]
    Pf_ld = [max(sol_laat.u[i][Model2.Pfind], 1.0) for i in eachindex(sol_laat.u)]
    B_ld  = [sol_laat.u[i][Model2.Benzind] for i in eachindex(sol_laat.u)]

    p2d_pop = plot(t_ld, [N_ld Lyt_ld Lys_ld Pf_ld],
        label=["Naïef N" "Lytisch L" "Lysogeen l" "Vrije fagen"],
        color=[:blue :red :green :darkred], lw=2,
        yscale=:log10, ylims=(1,:auto),
        xlabel="t [h]", ylabel="Cellen of fagen L⁻¹",
        legend=:topleft,
        bottom_margin=6Plots.mm, left_margin=10Plots.mm)
    vline!(p2d_pop, [t_infectie_laat], color=:black, lw=1.5,
        linestyle=:dash, label="t_infectie")

    p2d_benz = plot(t_ld, B_ld,
        label="Benzonase", color=:green, lw=2,
        xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
        legend=:topleft,
        bottom_margin=6Plots.mm, left_margin=10Plots.mm)

    fig2d = plot(p2d_pop, p2d_benz, layout=(1,2), size=(1000,430), margin=6Plots.mm)
    savefig(fig2d, "h2d_lage_beginbiomassa.png")
    println("  Figuur opgeslagen: h2d_lage_beginbiomassa.png")
end

# ============================================================
#  2e. Delta_t sweep: fagen VOOR maltosepuls toevoegen
#  Toont: (1) eindresultaat + (2) tijdstip van max Benzonase
# ============================================================
println("\n=== 2e: Fagen voor maltosepuls (delta_t sweep) ===")

function run_fagen_eerst(N0, t_inf, moi, delta_t_na, malt_puls)
    t_malt = t_inf + delta_t_na
    p = make_p2(N0, t_inf, moi)
    u0 = zeros(23)
    u0[Model2.Sind]  = [4.44, 0.0, 5.42, 0.0]
    u0[Model2.Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Model2.Nind]  = N0
    tspan = (0.0, p.duration)

    maltoseCondition(u,t,integrator)   = t == t_malt
    maltoseAffect!(integrator)         = integrator.u[Model2.Sind[2]] += malt_puls
    maltoseCallBack = DiscreteCallback(maltoseCondition, maltoseAffect!)
    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model2.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)
    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model2.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)
    domainCondition(u,t,integrator)    = any(x->x<0.0, u)
    domainAffect!(integrator)          = Model2.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(Model2.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([t_malt; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, maltoseCallBack,
                             infectionCallBack, fbaCallBack))
end

t_inf_2e           = 5.0
malt_2e            = 2.337
moi_2e             = 2.0
delta_t_na_waarden = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
max_benz_na        = Float64[]
t_max_benz_na      = Float64[]

for dt in delta_t_na_waarden
    sol_na = run_fagen_eerst(1e9, t_inf_2e, moi_2e, dt, malt_2e)
    if SciMLBase.successful_retcode(sol_na)
        B_na  = [sol_na.u[i][Model2.Benzind] for i in eachindex(sol_na.u)]
        i_max = argmax(B_na)
        push!(max_benz_na,   B_na[i_max])
        push!(t_max_benz_na, sol_na.t[i_max])
        println("  dt=$dt h na fagen: max_benz=$(round(B_na[i_max],sigdigits=3)) | t_max=$(round(sol_na.t[i_max],digits=2)) h")
    else
        push!(max_benz_na,   0.0)
        push!(t_max_benz_na, NaN)
    end
end

# Referentie: geen maltose
sol_geen_malt  = run_fagen_eerst(1e9, t_inf_2e, moi_2e, 999.0, 0.0)
B_geen         = [sol_geen_malt.u[i][Model2.Benzind] for i in eachindex(sol_geen_malt.u)]
benz_geen_malt = maximum(B_geen)
t_geen_max     = sol_geen_malt.t[argmax(B_geen)]

# Referentie: maltose voor fagen
benz_malt_voor = maximum(max_benz_glc)
t_malt_voor    = minimum(t_max_benz_glc)

dt_labels_e = string.(delta_t_na_waarden)

p2e_eind = plot(dt_labels_e, max_benz_na,
    marker=:circle, lw=2, color=:teal,
    xlabel="Δt maltose ná fagen [h]",
    ylabel="Max Benzonase [mmol L⁻¹]",
    label="Maltose na fagen",
    legend=:topright,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
hline!(p2e_eind, [benz_geen_malt],
    color=:red, lw=2, linestyle=:dash, label="Geen maltose")
hline!(p2e_eind, [benz_malt_voor],
    color=:black, lw=2, linestyle=:dot, label="Maltose vóór fagen (ref)")

p2e_tijdstip = plot(dt_labels_e, t_max_benz_na,
    marker=:square, lw=2, color=:darkorange,
    xlabel="Δt maltose ná fagen [h]",
    ylabel="Tijdstip max Benzonase [h]",
    label="Maltose na fagen",
    legend=:topright,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
hline!(p2e_tijdstip, [t_geen_max],
    color=:red, lw=2, linestyle=:dash, label="Geen maltose")
hline!(p2e_tijdstip, [t_malt_voor],
    color=:black, lw=2, linestyle=:dot, label="Maltose vóór fagen (ref)")

fig2e = plot(p2e_eind, p2e_tijdstip, layout=(1,2), size=(1000,430), margin=6Plots.mm)
savefig(fig2e, "h2e_fagen_voor_maltose.png")
println("  Figuur opgeslagen: h2e_fagen_voor_maltose.png")

# ============================================================
#  2f. Benzonase productie bij 4 MOI-waarden (Model 2)
# ============================================================
println("\n=== 2f: Benzonase productie bij verschillende MOI (Model 2) ===")

moi_benz    = [0.001, 0.1, 1.0, 5.0]
kleuren_2f  = [:steelblue :darkorange :green :darkred]
labels_2f   = ["MOI=0.001" "MOI=0.1" "MOI=1.0" "MOI=5.0"]

p2f_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
p2f_lys  = plot(xlabel="t [h]", ylabel="Lysogene cellen L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)

benz_eindwaarden = Float64[]
for (idx, moi) in enumerate(moi_benz)
    sol_f = Model2.run(make_p2(1e9, 2.0, moi))
    if SciMLBase.successful_retcode(sol_f)
        t_f   = sol_f.t
        B_f   = [sol_f.u[i][Model2.Benzind] for i in eachindex(sol_f.u)]
        Lys_f = max.([sol_f.u[i][Model2.lind] for i in eachindex(sol_f.u)], 1.0)
        push!(benz_eindwaarden, maximum(B_f))
        plot!(p2f_benz, t_f, B_f,   label=labels_2f[idx], color=kleuren_2f[idx], lw=2)
        plot!(p2f_lys,  t_f, Lys_f, label=labels_2f[idx], color=kleuren_2f[idx], lw=2)
    end
end

p2f_bar = bar(string.(moi_benz), benz_eindwaarden,
    xlabel="MOI [-]", ylabel="Max Benzonase [mmol L⁻¹]",
    color=[:steelblue :darkorange :green :darkred],
    legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)

fig2f = plot(p2f_lys, p2f_benz, p2f_bar, layout=(1,3), size=(1300,430), margin=6Plots.mm)
savefig(fig2f, "h2f_benzonase_moi_vergelijking.png")
println("  Figuur opgeslagen: h2f_benzonase_moi_vergelijking.png")

println("\n=== Hoofdstuk 2 voltooid ===")