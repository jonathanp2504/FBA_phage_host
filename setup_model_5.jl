# ============================================================
#  SETUP MODEL 5
#  Model 5 verschilt van Model 3 door een extra fbaModelLytic
#  voor lytische cellen (objectief: max faagproductie).
#  De Parameters struct heeft fbaModelLytic, phage_id en
#  q_phage_L als extra velden.
#
#  Gebruik:
#    include("setup_model_5.jl")
#    include("calibrate_refs.jl")
#    include("optimization_weighted_A.jl")
#    include("optimization_weighted_B.jl")
#    resultaat = run_optimization_B(p)
# ============================================================

using OrdinaryDiffEq
using DelayDiffEq
import OrdinaryDiffEqCore
if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end

include("./Bin_model_5/parameters.jl")
include("./Bin_model_5/FBA.jl")
include("./Bin_model_5/dFBA.jl")

using COBREXA, AbstractFBCModels
import SBMLFBCModels

model_path   = joinpath(@__DIR__, "iJO1366.xml")
naiveModel   = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
lyticModel   = addPhage!(loadFBAmodel(model_path), phage_stoich)

alpha_syn      = 2.0;   beta_deg = 0.5
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 3.75, 0.0, 4.0]
E_coli_cellDW  = 1.0e-12
infection_time = 2.0
K_mal          = 0.01
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
all_ex_ids     = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 6.0e-12
duration       = 20.0
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)

naiveFba   = buildFbaCache(naiveModel,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                           benz_id="R_BENZ_prod")
lyticFba   = buildLyticFbaCache(lyticModel, exchange_ids,
                                "R_BIOMASS_Ec_iJO1366_core_53p95M", "R_PHAGE_prod")

# Let op: fbaModelLytic, phage_id en q_phage_L als extra velden
p = Parameters(
    duration, 1e9,
    alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau,
    E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba, lysogenFba, lyticFba,  # drie FBA caches
    0.0, zeros(4),
    0.0, zeros(4),
    0.0,                # q_benz_l
    0.0,                # q_phage_L
    1e-10,              # k_attach
    10.0,               # k_dettach
    5.0,                # k_inject
    K_mal,
    infection_time,
    1e8,                # infection_dose
    "R_BENZ_prod",      # benz_id
    "R_PHAGE_prod",     # phage_id
    0.05,               # k_tox
    0.1,                # beta_benz
    0.0,                # q_benz
    mu_max_vector,
    e_max_vector,
    0.0015              # f_prod
)

println("Model 5 setup voltooid — p aangemaakt")
