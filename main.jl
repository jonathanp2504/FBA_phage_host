using Plots
using UnPack
using COBREXA
using AbstractFBCModels
using DelayDiffEq
using OrdinaryDiffEq
using SciMLBase
import SBMLFBCModels
include("./Bin/parameters.jl")   # parameters_model5.jl hernoemen
include("./Bin/FBA.jl")          # FBA_model5.jl hernoemen
include("./Bin/dFBA.jl")         # dFBA_model5.jl hernoemen
include("optimization_2.jl")

model_path = joinpath(@__DIR__, "iJO1366.xml")
@assert isfile(model_path) "iJO1366.xml does not exist at path: $model_path"

naiveModel   = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
lyticModel   = addPhage!(loadFBAmodel(model_path), phage_stoich)

alpha_syn      = 0.2*10
beta_deg       = 0.05*10
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 3.75, 0.0, 4.0]
E_coli_cellDW  = 1.0e-12
infection_time = 1.0
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e",
                  "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e",
                  "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e",
                  "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e",
                  "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 6.0e-12
duration       = 20.0
startingBiomass = 1e9
K_mal          = 0.01
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)

naiveFba   = buildFbaCache(naiveModel,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
lyticFba   = buildLyticFbaCache(lyticModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M", "R_PHAGE_prod")

p = Parameters(
    duration, startingBiomass,
    alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba, lysogenFba, lyticFba,
    0.0, zeros(4),
    0.0, zeros(4), 0.0,
    0.0,              # q_phage_L startwaarde
    1e-10, 10.0, 5.0, K_mal,
    infection_time, 1e8,
    "R_BENZ_prod",
    "R_PHAGE_prod",
    0.05, 0.1, 0.0,
    mu_max_vector, e_max_vector,
    0.0015
)

sol = run(p)
include("./plotting.jl")

println("--- Optimalisatie Model 5 gestart ---")
tspan   = (0.0, duration)
opt_sol = run_optimization(p, tspan)

best_t_inf   = opt_sol.best_t_inf
best_biomass = opt_sol.best_biomass

println("MOI (vast)         : 2.0")
println("Beste infectietijd : ", round(best_t_inf,   digits=2), " uur")
println("Beste beginbiomassa: ", round(best_biomass, digits=0), " cellen/L")
println("Max Benzonase      : ", round(opt_sol.best_benzonase, digits=6))

p_optimal = Parameters(
    p.duration, best_biomass,
    p.alpha_syn, p.beta_deg, p.K_s, p.V_max, p.p_pref,
    p.tau, p.E_coli_cellDW, p.MW, p.h_release,
    p.biomass_id, p.ex_ids, p.all_exchanges, p.essentials,
    p.fbaModelNaive, p.fbaModelLysogen, p.fbaModelLytic,
    p.mu_N, copy(p.q_N),
    p.mu_l, copy(p.q_l), p.q_benz_l,
    p.q_phage_L,
    p.k_attach, p.k_dettach, p.k_inject, p.K_mal,
    best_t_inf, 2.0 * best_biomass,
    p.benz_id, p.phage_id,
    p.k_tox, p.beta_benz, p.q_benz,
    copy(p.mu_max), copy(p.e_max), p.f_prod
)

sol_final = run(p_optimal)
plotAll(sol_final, p_optimal)
