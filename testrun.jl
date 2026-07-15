using Pkg
Pkg.activate(".")
include("Model2.jl")
using Plots, Statistics, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

gr()

const BENZ_REF_3 = 0.001
const PN_REF_3   = 2.0e10

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

tau_death = 1.0   # [h] -- fixed delay between production start and cell death (replaces k_tox)

naiveModel2   = Model2.loadFBAmodel(model_path)
all_ex_ids    = [id for id in keys(naiveModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2     = Model2.buildFbaCache(naiveModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2   = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                      "R_BIOMASS_Ec_iJO1366_core_53p95M")

function make_p2(N0, t_inf, moi, f_prod_val=0.0015)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", tau_death, 0.0, 0.0,
        mu_max_vector, e_max_vector, f_prod_val, 0.05)
end


p_mu = make_p2(1e9, 2.0, 1e-5)
sol_mu = Model2.run(p_mu)

X_tot = [Model2.getTotalBiomass(u, p_mu) for u in sol_mu.u]
plot(sol_mu.t, X_tot, yscale=:log10, ylims=(1, :auto), label="Totaal X",
     color=:black, lw=3, xlabel="t [h]", ylabel="cells [1/L]", title="Bacteriële Populatie")
plot!(sol_mu, idxs=[Model2.Nind, Model2.Dind, Model2.Lind, Model2.lind],
      label=["Naive (N)" "Deciding (D)" "Lytic (I)" "Lysogeen (L)"], alpha=0.7, legend=:bottomright)