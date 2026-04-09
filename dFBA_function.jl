
function dFBA_phage_system(du, u, h, p::Parameters, t)
    uDecision = h(p, t - 20/60) #
    uLysis = h(p, t - 60/60)
    #println("t = $t, mu = $(p.mu), X_tot = $(u[Nind])")
    # --- 1. EXTRACTIE ---
    S_subs = u[Sind]
    e_enz  = u[Eind]
    C_benz = u[Benzind]
    e_mal  = u[Eind[2]] # LamB proxy (Maltose enzym)

    # --- 2. CYBERNETICA & FBA ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    u_cyt = getU_cyt(R, p)
    v_cyt = getV_cyt(R)

    X_tot = getTotalBiomass(u, p)

    # Naive cells
    du[Nind] = p.mu * u[Nind] # growth
    du[Nind] -= p.k_inject * u[Paind] * u[Nind]/X_tot # infection --> N -> D

    # Deciding cells 
    du[Dind] = p.k_inject * u[Paind] * u[Nind]/X_tot # N -> D
    du[Dind] -= p.k_inject * uDecision[Paind] * uDecision[Nind]/getTotalBiomass(uDecision, p) # D -> L of l

    # Lytic cells
    du[Lind] = (1-getProbLys(u, p)) * p.k_inject * uDecision[Paind] * uDecision[Nind]/getTotalBiomass(uDecision, p) # D -> L
    du[Lind] -= (1-getProbLys(uDecision, p)) * p.k_inject * uLysis[Paind] * uLysis[Nind]/getTotalBiomass(uLysis, p) # L -> (lysis)

    # Lysogenic cells
    du[lind] = getProbLys(u, p) * p.k_inject * uDecision[Paind] * uDecision[Nind]/getTotalBiomass(uDecision, p) # D -> l
    du[lind] += p.mu * u[lind] # growth of lysogens

    # MOI 
    du[MOIind] = p.k_inject * u[Paind] / X_tot
    du[MOIind] -= p.k_inject * uDecision[Paind] / getTotalBiomass(uDecision, p)

    # Free phages
    du[Pfind] = p.b * (1-getProbLys(uDecision, p)) * p.k_inject * uLysis[Paind] * uLysis[Nind]/getTotalBiomass(uLysis, p)
    du[Pfind] -= p.k_attach * X_tot * u[Pfind] # attachment
    du[Pfind] += p.k_dettach * u[Paind] # dettachment

    # Attached phages
    du[Paind] = p.k_attach * X_tot * u[Pfind] # attachment
    du[Paind] -= p.k_dettach * u[Paind] # dettachment
    du[Paind] -= p.k_inject * u[Paind] # injection   
    
    # 4. DIFFERENTIAALVERGELIJKINGEN  
    for i in eachindex(Sind)
        sub_idx = Sind[i]
        enz_idx = Eind[i]

        # VERBETERDE LOGICA:
        # We laten de opname altijd toe, tenzij het substraat écht op is EN 
        # de FBA geen opname meer voorschrijft (p.q[i] >= 0).
        # Als er door lysis weer glucose bijkomt (S_subs[i] > 1e-7), 
        # zal de term p.q[i] * ... weer negatief worden en de glucose doen dalen.
        
        #if S_subs[i] < 1e-8 && p.q[i] >= 0
             #du[sub_idx] = 0.0
        #else
             #du[sub_idx] = p.q[i] * p.E_coli_cellDW * X_tot 
        #end
        # Voorkom dat du[sub_idx] de concentratie negatief maakt
        flux_val = p.q[i] * p.E_coli_cellDW * X_tot
        du[sub_idx] = max(0.0, flux_val)
        
        # Glucose vrijgave bij lysis (enkel voor index 1 = glucose)
        if i == 1
            du[sub_idx] += p.h_release * p.k_inject * uLysis[Paind] * uLysis[Nind]/getTotalBiomass(uLysis, p)
        end

        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - p.beta_deg* e_enz[i] + 0.001 # kleine basale expressie zodat er altijd een beetje enzym is (voor snelle adaptatie bij nieuwe substraat beschikbaar)
    end

    # --- 4. BENZONASE BALANS ---
    # Productie (via p.q_benz uit FBA) minus afbraak
    du[Benzind] = (p.q_benz * p.E_coli_cellDW * X_tot) #- (p.beta_benz * C_benz)
end

function getMonod(substrates::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return [max(0.0, substrates[i] / (substrates[i] + parameters.K_s[i])) for i in eachindex(substrates)]
end

function getRate(monod::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return parameters.V_max .* monod
end

# Berekent de cybernetische variabele u (synthese)
function getU_cyt(rate::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    weighted_R = parameters.p_pref .* rate
    denom = sum(weighted_R) + 1e-10
    return weighted_R ./ denom
end

# Berekent de cybernetische variabele v (activiteit)
function getV_cyt(rate::Vector{Float64})::Vector{Float64}
    denom = maximum(rate) + 1e-10
    return rate ./ denom
end

# Haalt veilig de fluxwaarden op voor de gevraagde IDs. 
# Als een ID niet in de oplossing zit, wordt 0.0 teruggegeven.
function getFluxes(sol, ids::Vector{String})::Vector{Float64}
    flux_vector = zeros(length(ids)) 
    for i in eachindex(ids)
        id = ids[i]
        if haskey(sol.fluxes, id)
            flux_vector[i] = sol.fluxes[id]
        else
            flux_vector[i] = 0.0
        end
    end
    return flux_vector
end