# ============================================================
#  RUN SCRIPT — OPTIMIZER MODEL 4
#  Optimizer A: vaste biomassa (1e9), optimaliseert MOI + t_inf
#  Optimizer B: vaste MOI (2.0), optimaliseert t_inf + biomassa
# ============================================================
include("Model4.jl")
include("optimizer_utils.jl")
include("run_optimalisatie_A_M4.jl")
include("run_optimalisatie_B_M4.jl")

using Plots, Printf, JuMP, HiGHS
using COBREXA, AbstractFBCModels
import SBMLFBCModels

model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 1.32
b0             = 50.0;  kb = 228.6
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

naiveModel4   = Model4.loadFBAmodel(model_path)
lysogenModel4 = Model4.addBenzonase!(Model4.loadFBAmodel(model_path), Model4.benz_stoich)
all_ex_ids    = [id for id in keys(naiveModel4.reactions) if startswith(id, "R_EX_")]
naiveFba4     = Model4.buildFbaCache(naiveModel4, exchange_ids,
                    "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba4   = Model4.buildFbaCache(lysogenModel4, exchange_ids,
                    "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")

p_initial = Model4.Parameters(duration, 1e9, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, b0, kb, 1e-12, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba4, lysogenFba4,
    0.0, zeros(4), 0.0, zeros(4), 0.0,
    7.92e-8, 6.48, 3.02, 0.01,
    2.0, 2.0*1e9,
    "R_BENZ_prod", 0.05, 0.1, 0.0,
    mu_max_vector, e_max_vector, 0.0015)

println("=== BENZ_REF berekening voor Model 4 ===")
println("\n" * "="^60)
println("OPTIMIZER A — MODEL 4")
println("="^60)
res_A = run_optimization_A(p_initial; w_benz=0.8, w_P=0.2)
println("\nResultaat A: MOI=$(round(res_A.best_moi,digits=4)) | ",
        "t_inf=$(round(res_A.best_t_inf,digits=2)) h | ",
        "score=$(round(res_A.best_score,digits=6))")

println("\n" * "="^60)
println("OPTIMIZER B — MODEL 4")
println("="^60)
res_B = run_optimization_B(p_initial; w_benz=0.8, w_P=0.2)
println("\nResultaat B: t_inf=$(round(res_B.best_t_inf,digits=2)) h | ",
        "N0=$(round(res_B.best_N0,sigdigits=3)) cellen/L | ",
        "score=$(round(res_B.best_score,digits=6))")
