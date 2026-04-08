function dFBA_phage_system(du, u, h, p::Parameters, t)
 
    # --- 1. EXTRACTIE ---
    S_subs = u[p.ind_subs]
    e_enz  = u[p.ind_e]
    I_cell, L_cell, P_phage = u[p.ind_I], u[p.ind_L], u[p.ind_P]
    C_benz = u[p.ind_Benz]
    e_mal  = u[p.ind_e[2]] # LamB proxy (Maltose enzym)

    # --- 2. CYBERNETICA & FBA ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    u_cyt = getU_cyt(R, p)
    v_cyt = getV_cyt(R)

    fbaUpdate!(u, p) # Berekent p.mu en p.q_benz
    X_tot = getTotalBiomass(u, p)

    # --- 3. ADSORPTIE LADDER LOGICA (Sequential Model) ---
    # Effectieve capaciteit o.b.v. LamB (e_mal)
    eff_cap = p.max_phages * (e_mal / (p.K_mal + e_mal + 1e-15))
    
    total_new_infections = 0.0
    total_ads_loss = 0.0 # We houden bij hoeveel fagen er totaal adsorberen (verlies voor P, want niet meer vrij in medium) 

    # We loopen door de ladder van S0 t/m S_max
    for i in 0:p.max_phages
        idx = p.ind_S_ladder[i+1] 
        S_i = u[idx]

        # A. TOXICITEIT BEREKENEN VOOR DEZE TREDE
        # S-cellen zijn 10x minder gevoelig dan L-cellen volgens jouw code
        toxic_death_S_i = (p.k_tox * 0.1) * C_benz * S_i

        # B. BINDING (S_i -> S_i+1)
        v_bind = 0.0 # snelheid waarmee cellen extra faag binden en naar de volgende trede gaan
        if i < p.max_phages
            binding_multiplier = max(0.0, (eff_cap - i) / p.max_phages)
            v_bind = p.k_on * P_phage * S_i * binding_multiplier
            
            du[idx] -= v_bind
            du[idx+1] += v_bind 
            total_ads_loss += v_bind
        end

        # C. LOSKOPPELING (S_i -> S_i-1)
        v_off = 0.0 # snelheid waarmee fagen loskomen van cellen en terug naar de vorige trede gaan (want nog niet geinfecteerd dus reversible)
        if i > 0
            v_off = p.k_off * i * S_i
            du[idx] -= v_off
            du[idx-1] += v_off 
            total_ads_loss -= v_off
        end

        # D. INJECTIE (S_i -> I)
        v_ins = 0.0
        if i > 0
            v_ins = p.k_ins * i * S_i
            du[idx] -= v_ins
            total_new_infections += v_ins
        end

        # E. BALANS VOOR S_i (Groei - Dood door Benz)
        # We voegen hier de mu * S_i en de toxic_death_S_i samen
        du[idx] += p.mu * S_i - toxic_death_S_i
    end

    # --- 4. PHAGE-HOST INTERACTIE ---
    # Vertraagde Lysis (gebaseerd op succesvolle injecties uit het verleden)
    #u_p = h(p, t - p.tau)
    # Sommatie van k_ins * i * S_i op tijdstip t-tau
    #total_ins_past = sum([p.k_ins * i * u_p[p.ind_S_ladder[i+1]] for i in 1:p.max_phages])
    #lysis_term = t > p.tau ? total_ins_past : 0.0
    # Efficiëntere berekening van de delay term
    # 1. Definieer de lysis_term met het juiste type (Float64 of Dual)
    lysis_term = 0.0
    if t > p.tau
        u_p = h(p, t - p.tau)
    # Loop handmatig door de ladder voor snelheid en stabiliteit
        for i in 1:p.max_phages
        # Veiligheidscheck: zorg dat we niet buiten de u_p vector indexeren
        val_past = u_p[p.ind_S_ladder[i+1]]
        lysis_term += p.k_ins * i * max(0.0, val_past)
        end
    end

    # Lysogenie beslissing (gebaseerd op lokale MOI)
    phi = (total_new_infections / (X_tot + 1e-10)) * 0.1 # Jouw MOI factor
    prob_lys = getProbLys(phi)

    # Toxiciteit voor L-cellen
    toxic_death_L = p.k_tox * C_benz * L_cell

    # Balansen voor I, L en P
    du[p.ind_I] = (1 - prob_lys) * total_new_infections - lysis_term
    du[p.ind_L] = p.mu * L_cell + (prob_lys * total_new_infections) - toxic_death_L
    du[p.ind_P] = (p.b * lysis_term) - total_ads_loss
    # 4. DIFFERENTIAALVERGELIJKINGEN  
    for i in 1:length(p.ind_subs)
        sub_idx = p.ind_subs[i]
        enz_idx = p.ind_e[i]

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
        if S_subs[i] <= 0.0 && flux_val > 0.0
            du[sub_idx] = 0.0
        else
            du[sub_idx] = flux_val
        end
        
        # Glucose vrijgave bij lysis (enkel voor index 1 = glucose)
        if i == 1
            du[sub_idx] += p.h_release * lysis_term
        end

        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - p.beta_deg* e_enz[i] + 0.001 # kleine basale expressie zodat er altijd een beetje enzym is (voor snelle adaptatie bij nieuwe substraat beschikbaar)
    end

    # --- 4. BENZONASE BALANS ---
    # Productie (via p.q_benz uit FBA) minus afbraak
    du[p.ind_Benz] = (p.q_benz * p.E_coli_cellDW * X_tot) #- (p.beta_benz * C_benz)
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