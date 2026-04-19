
using Plots
#Pkg.activate(".")
#Pkg.add("UnPack")
# Verwijder de Pkg.add regels uit je script als ze al geïnstalleerd zijn
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
naiveModel = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)

#model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].upper_bound = 2.0 # [1/h] moet bovenlimiet op staan anders te hoog
# 2. Cybernetische Parameters 
alpha_syn = 0.2*10   # Snelheid van enzymsynthese waardes van deze 2 nog opzoeken (stonden te laag)
beta_deg = 0.05*10   # Snelheid van enzymdegradatie
K_s = [0.0278, 0.0146, 0.0543, 0.0833]      # Affiniteit uit Tabel 7 (mmol/L)
tau = 1.0          # Latente periode (uren)
b = 170.0          # Burst size
alfa_ads = 1e-10  # Adsorptieconstante (L/gDW/h) lager dan bij Luan want eigenlijk meer deze waarde voor lambda bij te hoge alfa numerieke problemen), ik wil 2-step want moet voor hoge MOI zoals bij ons uiteindelijk gaat zijn
p_pref = [0.8925, 0.08925,  0.008925,  0.008925] # Voorkeurshiërarchie uit thesis
V_max = [12.7, 3.75, 0.0, 4.0]           # Opnamesnelheden (mmol/hr^-1)
E_coli_cellDW = 1.0e-12 # gDW per cel
infection_time = 1.0 # Tijdstip van infectie (uren)
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e", "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e", "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e", "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e", "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
naiveFba = buildFbaCache(naiveModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release = 6.0e-12 # glucose vrijgegeven bij lysis (mmol/burst) waarde gebaseerd op Luan table 7
duration = 20.0
startingBiomass = 1e9 # Start met 1 miljoen cellen per liter (1e6 cells/L)
n_hill = 2.0
# K_mu = 0.2
f_res = 0.05
K_mal = 0.01
mu_max_glc = 1.33
mu_max_mal = 1.26
mu_max_gly = 1.10
mu_max_ac = 0.29
mu_max_vector = [1.33, 1.26, 1.10, 0.29] 
# 2. Bereken e_max vector (Steady-state: synthese / (degradatie + groei))
# Formule: (alpha_syn + delta) / (beta_deg + mu_max)
e_max_vector = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)# glc, mal, glyc, ac
# --- 2. De Struct aanmaken ---
# Zorg dat de volgorde in je Parameters-file exact matcht met deze aanroep:
p = Parameters(
    duration, startingBiomass,
    # Cybernetica
    alpha_syn, beta_deg, K_s, V_max, p_pref, 
    
    # Phage-Host
    tau, b, E_coli_cellDW, MW_values, h_release, 

    # Model IDs
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids, naiveFba, lysogenFba,
    0.0,               # mu_N (startwaarde groei naïef)
    zeros(4),          # q_N  (startwaarde opname naïef)
    0.0,               # mu_l (startwaarde groei lysogeen)
    zeros(4),          # q_l  (startwaarde opname lysogeen)
    0.0,               # q_benz_l (startwaarde productie benzonase)

    # (Sequential Model)
    1e-10,          # k_on (jouw alfa_ads)
    10.0,           # k_off (1/h) -> Faag valt relatief snel los
    5.0,            # k_ins (1/h) -> Gemiddeld 12 min voor injectie
    K_mal,          # 0.01 (Verzadiging van LamB)

    infection_time, 
    1e8,

    "R_BENZ_prod", 
    0.05,           # k_tox
    0.1,            # beta_benz
    0.0,            # q_benz start op 0    
    mu_max_vector,  # DE NIEUWE VECTOR
    e_max_vector,    # DE NIEUWE VECTOR
    0.0015            #f_prod (% van de totale import wordt opgeofferd aan benzonase productie)
)


#sol = run(p)

# 6. UITGEBREID PLOTTEN
#include("./plotting.jl")
#plotAll(sol, p)
# --- 2. START OPTIMALISATIE ---
#println("--- Optimalisatie proces gestart ---")
#tspan = (0.0, duration)

# Roep de functie aan uit optimizer.jl
# Deze functie moet 'Optimization.solve' aanroepen en de beste x teruggeven
#opt_sol = run_optimization(p, tspan)

# Haal de winnende waarden op
#best_moi = opt_sol.u[1]
#best_t_inf = opt_sol.u[2]

#println("--- Optimalisatie voltooid ---")
#println("Beste MOI gevonden: ", round(best_moi, digits=4))
#println("Beste Infectietijd gevonden: ", round(best_t_inf, digits=2), " uur")
#println("Maximale Benzonase opbrengst: ", round(-opt_sol.objective, digits=6))

# --- 3. FINALE RUN MET OPTIMALE WAARDEN ---
# We overschrijven de initiële parameters met de resultaten van de optimizer
#p_optimal = remake(p, 
    #MOI = best_moi, 
    #infection_time = best_t_inf
#)

# Run de simulatie één keer met de beste instellingen
#sol_final = run(p_optimal)

# --- 4. PLOTTEN ---
# Nu plotten we de uitkomst van de geoptimaliseerde run
#plotAll(sol_final, p_optimal)

sol = run(p)

include("./plotting.jl")

# --- OPTIMALISATIE ---
println("--- Optimalisatie proces gestart ---")
tspan = (0.0, duration)

opt_sol = run_optimization(p, tspan)

best_t_inf    = opt_sol.best_t_inf
best_biomass  = opt_sol.best_biomass

println("MOI (vast)         : 2.0")
println("Beste infectietijd : ", round(best_t_inf, digits=2), " uur")
println("Beste beginbiomassa: ", round(best_biomass, digits=0), " cellen/L")
println("Max Benzonase      : ", round(opt_sol.best_benzonase, digits=6))

p_optimal = Parameters(
    p.duration,
    best_biomass,        # <-- variabel
    p.alpha_syn, p.beta_deg, p.K_s, p.V_max, p.p_pref,
    p.tau, p.b, p.E_coli_cellDW, p.MW, p.h_release,
    p.biomass_id, p.ex_ids, p.all_exchanges, p.essentials,
    p.fbaModelNaive, p.fbaModelLysogen,
    p.mu_N, copy(p.q_N), p.mu_l, copy(p.q_l), p.q_benz_l,
    p.k_attach, p.k_dettach, p.k_inject, p.K_mal,
    best_t_inf,
    2.0 * best_biomass,  # infectiedosis = MOI * N0
    p.benz_id, p.k_tox, p.beta_benz, p.q_benz,
    copy(p.mu_max), copy(p.e_max), p.f_prod
)

sol_final = run(p_optimal)
plotAll(sol_final, p_optimal)