function dFBA_phage_system(du, u, h, p::Parameters, t)
    println(t)
 
    # --- EXTRACTIE ---
    S_subs = u[p.ind_subs]
    e_enz  = u[p.ind_e]
    S_cell, I_cell, L_cell, P_phage = u[p.ind_S], u[p.ind_I], u[p.ind_L], u[p.ind_P]

    # --- 1. CYBERNETICA ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    u_cyt = getU_cyt(R, p)
    v_cyt = getV_cyt(R) # <--- De rode lijn zou nu moeten verdwijnen

   
    phi = getPhi(P_phage, p)
    prob_lys = getProbLys(phi)
    
    X_tot                 = getTotalBiomass(u, p)
    totaal_adsorptie_flux = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie_flux  = getNewInfectionFlux(u, p)
    u_p = h(p, t - p.tau)
    lysis_term = p.alfa_ads * u_p[p.ind_S] * u_p[p.ind_P]
    # 4. DIFFERENTIAALVERGELIJKINGEN  
    for i in 1:length(p.ind_subs)
        sub_idx = p.ind_subs[i] # Dit is index 1, 2, 3 of 4
        enz_idx = p.ind_e[i]    # Dit is index 5, 6, 7 of 8

        # Substraat verandering (mmol/L/h)
        du[sub_idx] = p.q[i] * p.E_coli_cellDW * X_tot 
        
        # Glucose release (alleen bij het eerste substraat)
        if i == 1
            du[sub_idx] += p.h_release * lysis_term
        end

        # Enzymdynamiek (gebruik de enzym_idx uit je parameters)
        # f[i], u_cyt[i] en e_enz[i] hebben index 1 t/m 4
        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - p.beta_deg * e_enz[i]
    end

    # Faag-Host interactie
    du[p.ind_S] = p.mu * S_cell - nieuwe_infectie_flux
    du[p.ind_I] = (1 - prob_lys) * nieuwe_infectie_flux - lysis_term
    du[p.ind_L] = p.mu * L_cell + (prob_lys * nieuwe_infectie_flux)
    du[p.ind_P] = (p.b * lysis_term) - totaal_adsorptie_flux
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