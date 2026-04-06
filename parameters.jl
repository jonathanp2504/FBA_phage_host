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
    alfa_ads::Float64
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

    # Indices voor state vector u
    ind_subs::Vector{Int64}
    ind_e::Vector{Int64}
    ind_S::Int64
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
end

# Berekent de lokale MOI (phi): de verhouding tussen infectiesnelheid en groeisnelheid
function getPhi(u, p::Parameters)
    alfa_eff = getAlfa_eff(u, p)
    P_phage = u[p.ind_P]
    return (alfa_eff * P_phage) * 0.1
end

# Berekent de kans op lysogenie op basis van phi
function getProbLys(phi::Float64)
    # Poisson n >= 2: 1 - P(0) - P(1)
    # P(0) = e^-phi, P(1) = phi * e^-phi
    return 1 - exp(-phi) - (phi * exp(-phi))
end

function getTotalBiomass(u, p::Parameters)
    return u[p.ind_S] + u[p.ind_I] + u[p.ind_L]
end

function getPhages(u, p::Parameters)::Float64 
    return u[p.ind_P]
end

function getTotalAdsorptionFlux(u, p::Parameters)
    X_tot = getTotalBiomass(u, p)
    alfa_eff = getAlfa_eff(u, p)
    return alfa_eff * X_tot * u[p.ind_P]
end

function getNewInfectionFlux(u, p::Parameters)
    alfa_eff = getAlfa_eff(u, p)
    return alfa_eff * u[p.ind_S] * u[p.ind_P]
end

function getAlfa_eff(u, p::Parameters)
    e_mal = u[p.ind_e[2]] # Stel dat maltose het 2e substraat is
    # e_mal loopt van 0 tot 1 (relatief enzymniveau)
    # We koppelen dit aan de maximale alfa_ads
    base_leak = 0.001 # 0.1% basale expressie
    return p.alfa_ads * (base_leak + (1 - base_leak) * e_mal)
end

function fbaUpdate!(u, p::Parameters)
    S_subs = u[p.ind_subs]
    e_enz  = u[p.ind_e]

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
                v_max_available = conc / (p.E_coli_cellDW * getTotalBiomass(u, p) + 1e-10)
                
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