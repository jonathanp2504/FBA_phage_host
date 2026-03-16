function dFBA_phage_system(du, u, h, p::Parameters, t)
    S_subs = u[p.ind_subs]
    e_enz  = u[p.ind_e]
    S_cell, I_cell, L_cell, P_phage = u[p.ind_S], u[p.ind_I], u[p.ind_L], u[p.ind_P]
    X_tot = getTotalBiomass(u, p)

    # 1. CYBERNETICA (f, u en v)
    f = getMonod(S_subs, p)
    R = V_max .* f
    
    # u_cyt bepaalt de synthese van enzymen (e)
    denom = sum(p_pref .* R) + 1e-10
    u_cyt = [(p_pref[i] * R[i]) / denom for i in 1:3]
    
    # v_cyt moduleert de werkelijke opname (v)
    v_cyt = [R[i] / (maximum(R) + 1e-10) for i in 1:3]

    # 2. LYSIS TERMIJN 
    u_p = h(p, t - tau)
    # Cellen die nu barsten zijn tau geleden geïnfecteerd
    lysis_term = p.alfa_ads * u_p[p.ind_S] * u_p[p.ind_P]

    # 3. FBA KOPPELING
    # for id in keys(model.reactions); if startswith(id, "R_EX_"); model.reactions[id].lower_bound = 0.0; end; end
    
    # De bounds worden beperkt door het huidige enzym-niveau e_enz
    model.reactions["R_EX_glc__D_e"].lower_bound = -R[1] * e_enz[1] * v_cyt[1]
    model.reactions["R_EX_mal__L_e"].lower_bound = -R[2] * e_enz[2] * v_cyt[2]
    model.reactions["R_EX_ac_e"].lower_bound      = -R[3] * e_enz[3] * v_cyt[3]
    
    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    
    if isnothing(sol) || isnan(sol.fluxes["R_BIOMASS_Ecoli_core_w_GAM"])
        mu, q = 0.0, zeros(3)
    else
        mu = sol.fluxes["R_BIOMASS_Ecoli_core_w_GAM"]
        q = [sol.fluxes["R_EX_glc__D_e"], sol.fluxes["R_EX_mal__L_e"], sol.fluxes["R_EX_ac_e"]]
    end

    # 4. DIFFERENTIAALVERGELIJKINGEN  
    for i in p.ind_subs
        du[i] = q[i] * p.E_coli_cellDW * X_tot 
        # Enzymdynamiek: Synthese (alfa * f * u) - Degradatie (beta * e)
        du[i+3] = p.alpha_syn * f[i] * u_cyt[i] - p.beta_deg * e_enz[i]
    end

    # Faag-Host interactie
    du[p.ind_S] = mu * S_cell - (p.alfa_ads * S_cell * P_phage)
    du[p.ind_I] = (p.alfa_ads * S_cell * P_phage) - lysis_term
    du[p.ind_L] = mu * L_cell
    du[p.ind_P] = (p.b * lysis_term) - (p.alfa_ads * S_cell * P_phage)
    println(t)
end

function getMonod(substrates::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return [max(0.0, substrates[i] / (substrates[i] + parameters.K_s[i])) for i in eachindex(substrates)]
end
