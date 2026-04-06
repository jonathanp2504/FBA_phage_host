function dFBA_phage_system(du, u, h, p::Parameters, t)
 
    # --- EXTRACTIE ---
    S_subs = u[p.ind_subs]
    e_enz  = u[p.ind_e]
    S_cell, I_cell, L_cell, P_phage = u[p.ind_S], u[p.ind_I], u[p.ind_L], u[p.ind_P]
    C_benz = u[p.ind_Benz] # De 13e variabele in u is de benzonase concentratie
    # --- 1. CYBERNETICA ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    u_cyt = getU_cyt(R, p)
    v_cyt = getV_cyt(R) # <--- De rode lijn zou nu moeten verdwijnen

   
    phi = getPhi(u, p)
    prob_lys = getProbLys(phi)
    
    X_tot                 = getTotalBiomass(u, p)
    totaal_adsorptie_flux = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie_flux  = getNewInfectionFlux(u, p)
    u_p = h(p, t - p.tau)
    lysis_term = t > p.tau ? getNewInfectionFlux(u_p, p) : 0.0
    mu_safe = max(0.0, p.mu) # Zorg dat groeisnelheid niet negatief wordt, met toevoeging van mu_safe crasht simulatie
    # --- 2. TOXICITEIT BEREKENEN ---
    # We berekenen de sterfte door Benzonase
    toxic_death_L = p.k_tox * C_benz * L_cell
    toxic_death_S = (p.k_tox * 0.1) * C_benz * S_cell # S is minder gevoelig
    # 4. DIFFERENTIAALVERGELIJKINGEN  
    for i in 1:length(p.ind_subs)
        sub_idx = p.ind_subs[i]
        enz_idx = p.ind_e[i]

        # VERBETERDE LOGICA:
        # We laten de opname altijd toe, tenzij het substraat écht op is EN 
        # de FBA geen opname meer voorschrijft (p.q[i] >= 0).
        # Als er door lysis weer glucose bijkomt (S_subs[i] > 1e-7), 
        # zal de term p.q[i] * ... weer negatief worden en de glucose doen dalen.
        
        if S_subs[i] < 1e-8 && p.q[i] >= 0
             du[sub_idx] = 0.0
        else
             du[sub_idx] = p.q[i] * p.E_coli_cellDW * X_tot 
        end
        
        # Glucose vrijgave bij lysis (enkel voor index 1 = glucose)
        if i == 1
            du[sub_idx] += p.h_release * lysis_term
        end

        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - p.beta_deg* e_enz[i] + 0.001 # kleine basale expressie zodat er altijd een beetje enzym is (voor snelle adaptatie bij nieuwe substraat beschikbaar)
    end

    # Faag-Host interactie
    du[p.ind_S] = p.mu * S_cell - nieuwe_infectie_flux - toxic_death_S
    du[p.ind_I] = (1 - prob_lys) * nieuwe_infectie_flux - lysis_term
    du[p.ind_L] = p.mu * L_cell + (prob_lys * nieuwe_infectie_flux) - toxic_death_L # raar dat er geen overshoot is van total X wat vreemd is, waarschijnlijk iets te maken met getphi ---> bekijk voor alternatief (misschien op andere manier doen, zonder growth speed initieel)
    du[p.ind_P] = (p.b * lysis_term) - totaal_adsorptie_flux 
    # --- 4. BENZONASE BALANS ---
    # Productie (via p.q_benz uit FBA) minus afbraak
    du[p.ind_Benz] = (p.q_benz * p.E_coli_cellDW * X_tot) - (p.beta_benz * C_benz)
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