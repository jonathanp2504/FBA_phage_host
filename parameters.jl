mutable struct Parameters 
    # Cybernetica & Kinetiek
    alpha_syn::Float64 
    beta_deg::Float64
    K_s::Vector{Float64}
    V_max::Vector{Float64}
    p_pref::Vector{Float64}
    
    # Phage-Host Parameters
    tau::Float64
    b::Float64
    E_coli_cellDW::Float64
    MW::Vector{Float64}
    h_release::Float64

    # Model IDs (iJO1366 specifiek)
    biomass_id::String
    ex_ids::Vector{String}
    all_exchanges::Vector{String}
    essentials::Vector{String}
    fbaModel
    mu::Float64
    q::Vector{Float64}

    # --- NIEUWE ADSORPTIE PARAMETERS (LADDER) ---
    max_phages::Int64    # Hoeveel fagen kunnen er op 1 cel? (bijv. 10)
    k_on::Float64        # Binding rate (L / gDW / h) - voorheen alfa_ads
    k_off::Float64       # Loskoppelingsrate (1/h)
    k_ins::Float64       # Injectierate (1/h) - de stap van S_i naar I
    K_mal::Float64       # Affiniteit voor LamB (verzadiging van de ladder)

    # Indices voor state vector u
    ind_subs::Vector{Int64}
    ind_e::Vector{Int64}
    ind_S_ladder::UnitRange{Int} # S0, S1, ..., S10
    ind_I::Int64
    ind_L::Int64
    ind_P::Int64
    infection_time::Float64

    # benzonase parameters
    ind_Benz::Int      # De index in vector u
    benz_id::String    # Het ID van de reactie in de FBA
    k_tox::Float64     # Toxiciteitscoëfficiënt
    beta_benz::Float64 # Afbraaksnelheid enzym
    q_benz::Float64    # Opslag voor de berekende flux

    # Referentie mu_max waarden uit Luan
    mu_max_glc::Float64   # 1.33
    mu_max_mal::Float64   # 1.26
    mu_max_gly::Float64   # 1.10
    mu_max_ac::Float64    # 0.29
end


# 2. Receptor-factor (Storms: Michaelis-Menten verzadiging van LamB)
#function f_receptor(u, p::Parameters)
    # We halen het relatieve niveau van maltose-enzymen op
    #e_mal = u[p.ind_e[2]] 
    #return e_mal / (p.K_mal + e_mal + 1e-15)
#end


## --- POPULATIE DYNAMICA ---

function getPhi(u, p::Parameters)
    # Totaal aantal fagen gebonden aan cellen
    fagen_op_cellen = sum([i * u[p.ind_S_ladder[i+1]] for i in 1:p.max_phages])
    X_tot = getTotalBiomass(u, p)
    return fagen_op_cellen / (X_tot + 1e-15)
end

function getProbLys(phi::Float64)
    return 1 - exp(-phi) - (phi * exp(-phi))
end

function getTotalBiomass(u, p::Parameters)
    # Som van de ladder (S0 t/m S10) + geïnfecteerd + lysogeen
    return sum(u[p.ind_S_ladder]) + u[p.ind_I] + u[p.ind_L]
end

#function getTotalAdsorptionFlux(u, mu_werkelijk, u_cyt, p::Parameters)
    #X_tot = getTotalBiomass(u, p)
    #alpha_nu = get_alpha_eff(u, mu_werkelijk, u_cyt, p)
    #return alpha_nu * X_tot * u[p.ind_P]
#end

#function getNewInfectionFlux(u, mu_werkelijk, u_cyt, p::Parameters)
    #alpha_nu = get_alpha_eff(u, mu_werkelijk, u_cyt, p)
    #return alpha_nu * u[p.ind_S] * u[p.ind_P]
#end

#function getAlfa_eff(u, p::Parameters)
    #e_mal = u[p.ind_e[2]] # Stel dat maltose het 2e substraat is
    #e_mal loopt van 0 tot 1 (relatief enzymniveau)
    #We koppelen dit aan de maximale alfa_ads
    #base_leak = 0.001 # 0.1% basale expressie
    #return p.alfa_ads * (base_leak + (1 - base_leak) * e_mal)
#end

function fbaUpdate!(u, p::Parameters)
    S_subs = u[p.ind_subs]
    e_enz  = u[p.ind_e]
    println("t = ???, mu = $(p.mu), X_tot = $(getTotalBiomass(u, p))")
    # --- 1. CYBERNETISCHE BEREKENINGEN ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    v_cyt = getV_cyt(R)

    # --- 2. CONSTRAINTS UPDATEN VOOR ELK SUBSTRAAT ---
    for i in 1:length(p.ind_subs)
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
    L_cell = u[p.ind_L]
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
        p.q_benz = sol.fluxes[p.benz_id]
    else 
        # Als de cel niet kan groeien (bijv. geen suikers), mu op 0
        p.mu = 0.0
        p.q .= 0.0
        p.q_benz = 0.0
    end
end

function getPhages(u, p::Parameters)
    return u[p.ind_P]
end

function getTotalPhages(u, p::Parameters)
    free_phages = u[p.ind_P]
    bound_phages = sum([i * u[p.ind_S_ladder[i+1]] for i in 1:p.max_phages])
    return free_phages + bound_phages
end