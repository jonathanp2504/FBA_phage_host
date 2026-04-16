include("./Bin/dFBA.jl")


model_path = joinpath(@__DIR__, "iJO1366.xml")
t_inf_as = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]           # 6 stappen i.p.v. 12
moi_as   = [0.0001, 0.001, 0.01, 0.1, 1.0]      # We laten 100.0 weg (optioneel, voor extra snelheid)
heatmap_results = zeros(length(moi_as), length(t_inf_as))


println("Starten van de parameter sweep...")

# --- 1. VOORBEREIDING ---
duration = 15.0
naiveModel = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
startingBiomass = 1e6 # Start met 1 miljoen cellen per liter (1e6 cells/L)
alpha_syn = 0.2*10   # Snelheid van enzymsynthese waardes van deze 2 nog opzoeken (stonden te laag)
beta_deg = 0.05*10   # Snelheid van enzymdegradatie
K_s = [0.0278, 0.0146, 0.0543, 0.0833]      # Affiniteit uit Tabel 7 (mmol/L)
tau = 1.0          # Latente periode (uren)
b = 170.0          # Burst size
alfa_ads = 1e-10  # Adsorptieconstante (L/gDW/h) lager dan bij Luan want eigenlijk meer deze waarde voor lambda bij te hoge alfa numerieke problemen), ik wil 2-step want moet voor hoge MOI zoals bij ons uiteindelijk gaat zijn
p_pref = [0.8925, 0.08925,  0.008925,  0.008925] # Voorkeurshiërarchie uit thesis
V_max = [12.7, 3.75, 0.0, 4.0]           # Opnamesnelheden (mmol/hr^-1)
E_coli_cellDW = 1.0e-12 # gDW per cel
infection_time = 3.0 # Tijdstip van infectie (uren)
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e", "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e", "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e", "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e", "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
naiveFba = buildFbaCache(naiveModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release = 0.0 #6.0e-12 # glucose vrijgegeven bij lysis (mmol/burst) waarde gebaseerd op Luan table 7
duration = 20.0
n_hill = 2.0
# K_mu = 0.2
f_res = 0.05
K_mal = 0.01
mu_max_glc = 1.33
mu_max_mal = 1.26
mu_max_gly = 1.10
mu_max_ac = 0.29
mu_max_vector = [mu_max_glc, mu_max_mal, mu_max_gly, mu_max_ac] 
# 2. Bereken e_max vector (Steady-state: synthese / (degradatie + groei))
# Formule: (alpha_syn + delta) / (beta_deg + mu_max)
e_max_vector = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)# glc, mal, glyc, ac
# Zorg dat alle functies (fbaUpdate!, dFBA_phage_system) en indices geladen zijn.
# Gebruik de assen zoals je ze had gedefinieerd

for (i, moi_val) in enumerate(moi_as)
    for (j, t_inf_val) in enumerate(t_inf_as)
        println("\n[RUN] MOI: $moi_val | t_inf: $t_inf_val")
        run_label = "MOI=$(moi_val), t_inf=$(t_inf_val)"

        # 1. Maak een vers model voor deze specifieke iteratie
        naive_model_in_loop = loadFBAmodel(model_path)
        lysogen_model_in_loop = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
        naive_cache_in_loop = buildFbaCache(naive_model_in_loop, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
        lysogen_cache_in_loop = buildFbaCache(lysogen_model_in_loop, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")

        # 2. Bouw de parameters op met dit schone model
        p_loop = Parameters(
            duration, startingBiomass,
            alpha_syn, beta_deg, K_s, V_max, p_pref, 
            tau, b, E_coli_cellDW, MW_values, h_release, 
            "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
            naive_cache_in_loop, lysogen_cache_in_loop,
            0.0, zeros(4), 0.0, zeros(4), 0.0, 
            1e-10, 10.0, 5.0, K_mal, 
            t_inf_val, moi_val * startingBiomass, "R_BENZ_prod", 0.05, 0.1, 0.0, 
            mu_max_vector, e_max_vector, 0.0015
        )

        sol_loop = run(p_loop)

        # 5. Resultaat verwerking met de verbeterde check
        if SciMLBase.successful_retcode(sol_loop)
            max_val = maximum(sol_loop[Benzind, :])
            heatmap_results[i, j] = max_val
            println("   [OK] Max Benzonase: $max_val")
        else
            println("   [FAIL] Solver stopte met code: $(sol_loop.retcode)")
            heatmap_results[i, j] = 0.0
        end
    end
end

# --- STAP 3: HEATMAP PLOTTEN ---
# Als de heatmap_results gevuld zijn, zal dit nu werken zonder GKS errors
using Plots.Measures # Nodig voor de mm marges


# Gebruik een iets grotere epsilon om ook de allerkleinste waarden te vangen
log_data = log10.(heatmap_results .+ 1e-15)

# Bepaal het bereik voor de schaal
# We willen dat de schaal diep genoeg gaat om de verschillen links te zien
bovengrens = maximum(log_data)
ondergrens = bovengrens - 5 # Toon 5 ordes van grootte (bijv. van 10^-8 naar 10^-13)

p_heat = heatmap(t_inf_as, 1:length(moi_as), log_data,
    yticks = (1:length(moi_as), string.(moi_as)),
    xlabel = "Tijdstip van infectie [h]",
    ylabel = "Initiële MOI (P/N)",
    title  = "Benzonase Productie (Log10)",
    colorbar_title = "\n log10(mmol/L)",
    color  = :viridis,
    
    # FIX VOOR DE SCHAAL:
    clims = (ondergrens, bovengrens), 
    
    # FIX VOOR DE MARGES:
    right_margin = 20Plots.mm, # Nog meer ruimte voor de labels
    left_margin  = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    
    # Maak de getallen op de colorbar mooier (bijv. -9.0 ipv -9.0000000001)
    colorbar_formatter = x -> sprintf("%.1f", x) 
)

display(p_heat)
