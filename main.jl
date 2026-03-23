using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels
using Pkg
#Pkg.add("DifferentialEquations")
#Pkg.add("DelayDiffEq")
using Plots
using DifferentialEquations


# 1. Model laden
model_path = joinpath(@__DIR__, "iJO1366.xml")
@assert isfile(model_path) "iJO1366.xml does not exist at path: $model_path"
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model(model_path))

# 2. Cybernetische Parameters 
alpha_syn = 0.2   # Snelheid van enzymsynthese
beta_deg = 0.05   # Snelheid van enzymdegradatie
K_s = [0.0278, 0.0146, 0.0833]      # Affiniteit uit Tabel 7 (mmol/L)
tau = 0.6          # Latente periode (uren)
b = 170.0          # Burst size
alfa_ads = 1e-7  # Adsorptieconstante (L/gDW/h)
p_pref = [0.8925, 0.08925,  0.008925,  0.008925] # Voorkeurshiërarchie uit thesis
V_max = [15.0, 13.0, 11.0, 4.0]           # Opnamesnelheden (mmol/hr^-1)
E_coli_cellDW = 2.8e-13 # gDW per cel
infection_time = 5.0 # Tijdstip van infectie (uren)
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e", "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e", "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e", "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e", "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release = 5.66e-12 # glucose vrijgegeven bij lysis (mmol/burst) waarde gebaseerd op Luan table 7
# --- 2. De Struct aanmaken ---
# Zorg dat de volgorde in je Parameters-file exact matcht met deze aanroep:
include("./parameters.jl")
p = Parameters(alpha_syn, beta_deg, K_s, tau, b, alfa_ads, p_pref, V_max, E_coli_cellDW,MW_values, h_release, "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids, 9, 10, 11, 12, 5:8, 1:4, infection_time)
include("./dFBA_function.jl")
# 5. INITIALISATIE
# [Glc, Mal, Glyc, Ac, e_glc, e_mal, e_Glyc, e_ac, S, I, L, P] (subs in mmol/l)
u0 = [4.44, 2.337, 5.42, 0.0, 0.95, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0]
tspan = (0.0, 15.0) 

infectionCondition(u, t, integrator) = t == p.infection_time 
infectionAffect!(integrator) = integrator.u[parameters.ind_P] = 1000.0
infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

prob = DDEProblem(dFBA_phage_system, u0, (p,t)->u0, tspan, parameters)
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-6, tstops=[5.0], callback=infectionCallBack )

# 6. UITGEBREID PLOTTEN
# Plot A: Substraatverloop
p1 = plot(sol, idxs=parameters.ind_subs, title="Substraten (mmol/L)", 
          label=["Glucose" "Maltose" "Acetaat"], lw=2, ylabel="Conc.")

# Plot B: Enzymatische Adaptatie 
p2 = plot(sol, idxs=parameters.ind_e, title="Enzym-niveaus (Cybernetische e)", 
          label=["e_Glc" "e_Mal" "e_Ac"], lw=2, ls=:dash, ylabel="Relatief niveau")


# Plot C: Host populatie
p3 = plot(sol.t, [getTotalBiomass(u, parameters) for u in sol.u], title="Bacteriële Populatie", 
          label="Totaal X", color=:black, lw=3, ylabel="gDW/L", ylims=:auto)
plot!(p3, sol, idxs=[parameters.ind_S, parameters.ind_I], label=["Vatbaar (S)" "Infect (I)"], alpha=0.7)


# Plot D: Faag Lambda Dynamiek
p4 = plot(sol, idxs=[parameters.ind_P], title="Fagen (P)", yscale=:log10, 
          color=:red, lw=2, label="Phage Lambda", ylabel="Log10(P)")


plot(p1, p2, p3, p4, layout=(2,2), size=(800,800), margin=5Plots.mm)