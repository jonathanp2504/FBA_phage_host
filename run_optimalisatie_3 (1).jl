# ============================================================
#  RUN SCRIPT — OPTIMIZER MODEL 3
#  Optimizer A: vaste biomassa (1e9), optimaliseert MOI + t_inf
#  Optimizer B: vaste MOI (2.0), optimaliseert t_inf + biomassa
# ============================================================
include("Model2.jl")
include("Model3.jl")
include("optimizer_utils.jl")
include("run_optimalisatie_A_M23.jl")


using Plots, Printf, JuMP, HiGHS
using COBREXA, AbstractFBCModels
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
h_release      = 1.71e-12;  duration = 40.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW  = 1.0e-12

naiveModel3   = Model3.loadFBAmodel(model_path)
lysogenModel3 = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
all_ex_ids    = [id for id in keys(naiveModel3.reactions) if startswith(id, "R_EX_")]
naiveFba3     = Model3.buildFbaCache(naiveModel3, exchange_ids,
                    "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3   = Model3.buildFbaCache(lysogenModel3, exchange_ids,
                    "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
function maak_params_M3(p_base, moi, t_inf)
    infectiedosis = moi * p_base.startingBiomass
    return Model3.Parameters(
        p_base.duration, p_base.startingBiomass,
        p_base.alpha_syn, p_base.beta_deg, p_base.K_s, p_base.V_max, p_base.p_pref,
        p_base.tau, p_base.b, p_base.E_coli_cellDW, p_base.MW, p_base.h_release,
        p_base.biomass_id, p_base.ex_ids, p_base.all_exchanges, p_base.essentials,
        p_base.fbaModelNaive, p_base.fbaModelLysogen,
        p_base.mu_N, copy(p_base.q_N), p_base.mu_l, copy(p_base.q_l), p_base.q_benz_l,
        p_base.k_attach, p_base.k_dettach, p_base.k_inject, p_base.K_mal,
        t_inf, infectiedosis,
        p_base.benz_id, p_base.k_tox, p_base.beta_benz, p_base.q_benz,
        copy(p_base.mu_max), copy(p_base.e_max), p_base.f_prod)
end

p_initial = Model3.Parameters(duration, 1e9, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, b, 1e-12, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba3, lysogenFba3,
    0.0, zeros(4), 0.0, zeros(4), 0.0,
    7.92e-8, 6.48, 3.02, 0.01,
    2.0, 2.0*1e9,
    "R_BENZ_prod", 0.35, 0.1, 0.0,
    mu_max_vector, e_max_vector, 0.0015)

# BENZ_REF: lysogenModel3 heeft Benzonase al toegevoegd
# biomassareactie wordt intern op 0 gezet door bereken_benz_ref_fba
println("=== BENZ_REF berekening voor Model 3 ===")

println("\n" * "="^60)
println("OPTIMIZER A — MODEL 3")
println("="^60)
res_A = run_optimization_A(p_initial; w_benz=11.0, w_P=0.112, w_t=0.20)
println("\nResultaat A: MOI=$(round(res_A.best_moi,digits=4)) | ",
        "t_inf=$(round(res_A.best_t_inf,digits=2)) h | ",
        "score=$(round(res_A.best_score,digits=6))")


println("\n" * "="^60)
println("PARETO GRID SCAN — MODEL 3")
println("="^60)
# DOOR dit:
moi_values   = [0.001, 0.005, 0.01, 0.05, 0.1, 0.3,
                0.5, 1.0, 2.0, 3.0, 5.0]
t_inf_values = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]
pareto_grid = pareto_grid_scan_fixed(
    p_initial, Model3.run,
    Model3.Benzind, Model3.Pfind, Model3.Nind,
    moi_values, t_inf_values,
    maak_params_M3; figbase="pareto_M3_fixed")