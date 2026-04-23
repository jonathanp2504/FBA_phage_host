using Plots
using UnPack
using COBREXA
using AbstractFBCModels
using DelayDiffEq
using OrdinaryDiffEq
using SciMLBase
import SBMLFBCModels
include("./Bin/parameters.jl")
include("./Bin/FBA.jl")
include("./Bin/dFBA.jl")
include("optimization_2.jl")

# 1. Model laden
model_path = joinpath(@__DIR__, "iJO1366.xml")
@assert isfile(model_path) "iJO1366.xml does not exist at path: $model_path"

model = loadFBAmodel(model_path)

# 2. Parameters
alpha_syn      = 0.2*10
beta_deg       = 0.05*10
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0
b              = 170.0
alfa_ads       = 1e-10
E_coli_cellDW  = 2.8e-13
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [15.0, 13.0, 11.0, 4.0]
infection_time = 5.0
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e",
                  "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e",
                  "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e",
                  "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e",
                  "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 5.66e-12
duration       = 15.0
startingBiomass = 1e6
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)

fbaCache = buildFbaCache(model, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")

p = Parameters(
    duration, startingBiomass,
    alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, b, alfa_ads, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    fbaCache,
    0.0, zeros(4),
    infection_time,
    1e3,            # infection_dose (fagen bij infectie)
    "R_BENZ_prod",
    0.05,           # k_tox
    0.1,            # beta_benz
    0.0,            # q_benz
    mu_max_vector,
    e_max_vector,
    0.2,            # f_prod: 20% van groei naar Benzonase
    5.0             # Y_benz: mmol Benzonase per gDW groei-verlies
)

sol = run(p)

# Plotten
p1 = plot(sol, idxs=collect(Sind),
          title="Substraten (mmol/L)",
          label=["Glucose" "Maltose" "Glycerol" "Acetaat"], lw=2,
          xlabel="t [h]", ylabel="Conc.")

p2 = plot(sol, idxs=collect(Eind),
          title="Enzym-niveaus",
          label=["e_Glc" "e_Mal" "e_Glyc" "e_Ac"], lw=2, ls=:dash,
          xlabel="t [h]", ylabel="Relatief niveau")

p3 = plot(sol.t, [getTotalBiomass(u, p) for u in sol.u],
          title="Bacteriële Populatie", yscale=:log10, ylims=(1, :auto),
          label="Totaal X", color=:black, lw=3, ylabel="cells [1/L]")
plot!(p3, sol, idxs=[Sind_S, Sind_I, Sind_L],
      label=["Vatbaar (S)" "Lytisch (I)" "Lysogeen (L)"], alpha=0.7)

p4 = plot(sol, idxs=Pfind,
          title="Fagen (P)", yscale=:log10, ylims=(1, :auto),
          color=:red, lw=2, label="Phage Lambda",
          xlabel="t [h]", ylabel="phages [1/L]")

p5 = plot(sol, idxs=Benzind,
          title="Benzonase Productie",
          color=:green, lw=2,
          xlabel="t [h]", ylabel="Conc. [mmol/L]",
          label="Benzonase")

plot(p1, p2, p3, p4, p5, layout=(3, 2), size=(900, 1000), margin=5Plots.mm)

# --- OPTIMALISATIE ---
println("--- Optimalisatie proces gestart ---")
tspan = (0.0, duration)

opt_sol = run_optimization(p, tspan)

best_t_inf   = opt_sol.best_t_inf
best_biomass = opt_sol.best_biomass

println("MOI (vast)         : 2.0")
println("Beste infectietijd : ", round(best_t_inf,   digits=2), " uur")
println("Beste beginbiomassa: ", round(best_biomass, digits=0), " cellen/L")
println("Max Benzonase      : ", round(opt_sol.best_benzonase, digits=6))

p_optimal = Parameters(
    p.duration,
    best_biomass,
    p.alpha_syn, p.beta_deg, p.K_s, p.V_max, p.p_pref,
    p.tau, p.b, p.alfa_ads, p.E_coli_cellDW, p.MW, p.h_release,
    p.biomass_id, p.ex_ids, p.all_exchanges, p.essentials,
    p.fbaModel,
    p.mu, copy(p.q),
    best_t_inf,
    2.0 * best_biomass,
    p.benz_id, p.k_tox, p.beta_benz, p.q_benz,
    copy(p.mu_max), copy(p.e_max),
    p.f_prod, p.Y_benz
)

sol_final = run(p_optimal)

p1f = plot(sol_final, idxs=collect(Sind),
           title="Substraten (mmol/L)",
           label=["Glucose" "Maltose" "Glycerol" "Acetaat"], lw=2)
p3f = plot(sol_final.t, [getTotalBiomass(u, p_optimal) for u in sol_final.u],
           title="Bacteriële Populatie", yscale=:log10, ylims=(1, :auto),
           label="Totaal X", color=:black, lw=3)
p5f = plot(sol_final, idxs=Benzind,
           title="Benzonase (optimaal)",
           color=:green, lw=2, label="Benzonase")
plot(p1f, p3f, p5f, layout=(2, 2), size=(900, 700), margin=5Plots.mm)