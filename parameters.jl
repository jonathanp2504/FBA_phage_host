struct Parameters 
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

    # Indices voor state vector u
    ind_subs::Vector{Int64}
    ind_e::Vector{Int64}
    ind_S::Int64
    ind_I::Int64
    ind_L::Int64
    ind_P::Int64
    infection_time::Float64
end

# Berekent de lokale MOI (phi): de verhouding tussen infectiesnelheid en groeisnelheid
function getPhi(P_phage, mu, p::Parameters)
    # phi = (alpha * P) / mu
    return (p.alfa_ads * P_phage) / (mu + 1e-10)
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

function getTotalAdsorptionFlux(u, p::Parameters)
    X_tot = getTotalBiomass(u, p)
    return p.alfa_ads * X_tot * u[p.ind_P]
end

function getNewInfectionFlux(u, p::Parameters)
    return p.alfa_ads * u[p.ind_S] * u[p.ind_P]
end