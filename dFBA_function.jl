
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
    du[Nind] = p.mu_N * u[Nind] # growth
    du[Nind] -= p.k_inject * u[Paind] * u[Nind]/X_tot # infection --> N -> D

    # Deciding cells 
    du[Dind] = p.k_inject * u[Paind] * u[Nind]/X_tot # N -> D
    du[Dind] -= p.k_inject * uDecision[Paind] * uDecision[Nind]/getTotalBiomass(uDecision, p) # D -> L of l

    # Lytic cells
    du[Lind] = (1-getProbLys(u, p)) * p.k_inject * uDecision[Paind] * uDecision[Nind]/getTotalBiomass(uDecision, p) # D -> L
    du[Lind] -= (1-getProbLys(uDecision, p)) * p.k_inject * uLysis[Paind] * uLysis[Nind]/getTotalBiomass(uLysis, p) # L -> (lysis)

    # Lysogenic cells
    du[lind] = getProbLys(u, p) * p.k_inject * uDecision[Paind] * uDecision[Nind]/getTotalBiomass(uDecision, p) # D -> l
    du[lind] += p.mu_l * u[lind] # growth of lysogens

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
        # 1. Bereken de flux zoals voorheen
        # Totale opname = (opname per N-cel * aantal N) + (opname per l-cel * aantal l)
        # We gebruiken hier u[Nind] en u[lind] in plaats van X_tot
        # We tellen alle cellen die nog metabool actief zijn mee:
        # 1. De gezonde cellen (N) + de cellen in transitie (D & L)
        actieve_vrije_biomassa = u[Nind] + u[Dind] + u[Lind]
        flux_N = p.q_N[i] * actieve_vrije_biomassa * p.E_coli_cellDW
    
        # 2. De lysogene cellen (l) met hun eigen (tragere) opnamesnelheid q_l
        flux_l = p.q_l[i] * u[lind] * p.E_coli_cellDW
    
        flux_totaal = flux_N + flux_l
        du[sub_idx] = min(0.0, flux_totaal)
        # 2. De zachte stop (Hill-factor)
        # K_safe zorgt dat de flux naar 0 gaat als de concentratie bijna 0 is.
        #K_safe = 1e-4 
        #soft_stop = u[sub_idx] / (u[sub_idx] + K_safe)

        # 3. Bereken de uiteindelijke verandering
        #du[sub_idx] = flux_val * soft_stop
        # Totaal aantal cellen voor de noemer
        total_cells = u[Nind] + u[Dind] + u[Lind] + u[lind]

        # Alleen N (naïef) en L (lysogeen) dragen bij aan de totale groei
        # D en I hebben mu = 0, dus die vallen weg uit de teller
        if total_cells > 1.0
            mu_avg = (p.mu_N * u[Nind] + p.mu_l * u[lind]) / total_cells
        else
            mu_avg = p.mu_N
        end
        # Glucose vrijgave bij lysis (enkel voor index 1 = glucose)
        if i == 1
            du[sub_idx] += p.h_release * p.k_inject * uLysis[Paind] * uLysis[Nind]/getTotalBiomass(uLysis, p)
        end

        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - (p.beta_deg + mu_avg)* e_enz[i] + 0.001 # kleine basale expressie zodat er altijd een beetje enzym is (voor snelle adaptatie bij nieuwe substraat beschikbaar)
    end

    # --- 4. BENZONASE BALANS ---
    # --- 5. BENZONASE PRODUCTIE ---
    # Alleen de lysogenen (u[lind]) dragen bij
    productie_benz = p.q_benz_l * p.E_coli_cellDW * u[lind]
    du[Benzind] = productie_benz - (p.beta_benz * u[Benzind])
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