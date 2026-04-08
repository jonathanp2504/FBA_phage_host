include("./parameters.jl")
using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels 
function fbaUpdate!(u, p::Parameters)
    S_subs = u[Sind]
    e_enz  = u[Eind]
   
    # --- 1. CYBERNETISCHE BEREKENINGEN ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    v_cyt = getV_cyt(R)

    # --- 2. CONSTRAINTS UPDATEN VOOR ELK SUBSTRAAT ---
    for i in eachindex(Sind)
        id = p.ex_ids[i]
        conc = S_subs[i]
        
        if haskey(p.fbaModel.reactions, id)
            # Check of er suiker aanwezig is (bijv. na lysis)
            if conc > 1e-7  
                # De theoretische opname op basis van kinetiek en enzymen
                # We gebruiken hier e_enz[i] (cybernetische e) en v_cyt
                v_max_kin = R[i] * (e_enz[i] + 1e-3) * v_cyt[i]
                
                # Numerieke veiligheid: nooit meer opnemen dan er fysiek is
                # (S_subs / (biomassa * dt)) - we houden het simpel:
                v_max_available =  conc / (p.E_coli_cellDW * getTotalBiomass(u, p) + 1e-10)
                
                # De lower_bound is negatief voor opname
                p.fbaModel.reactions[id].lower_bound = -min(v_max_kin, v_max_available)
            else
                # Als conc echt 0 is (of lager door afronding), mag er niks in
                p.fbaModel.reactions[id].lower_bound = 0.0
            end
        end
    end
    # --- NIEUW: BENZONASE PRODUCTIE FORCEEREN ---
    L_cell = u[lind]
    X_tot = getTotalBiomass(u, p)
    lysogen_fraction = L_cell / (X_tot + 1e-10)

    # We definiëren een productie-snelheid (mmol per gram bacterie per uur)
    # Hoe meer lysogenen, hoe hoger de collectieve 'vraag' in het FBA model
    k_prod = 2  # Tweak deze waarde voor meer/minder opbrengst
    
    if haskey(p.fbaModel.reactions, p.benz_id)
        # De lower_bound dwingt de FBA om Benzonase te maken
        p.fbaModel.reactions[p.benz_id].lower_bound = k_prod * lysogen_fraction
    end
    # --- 3. FBA OPLOSSEN ---
    # Gebruik een try-catch of check isnothing om crashes bij onoplosbaarheid te voorkomen
    sol = flux_balance_analysis(p.fbaModel, optimizer = HiGHS.Optimizer)

    if !isnothing(sol) && haskey(sol.fluxes, p.biomass_id)
        p.mu = sol.fluxes[p.biomass_id]
        p.q  = getFluxes(sol, p.ex_ids)
        #p.q_benz = sol.fluxes[Benzind]
    else 
        # Als de cel niet kan groeien (bijv. geen suikers), mu op 0
        p.mu = 0.0
        p.q .= 0.0
        p.q_benz = 0.0
    end
end