# indices in state vector u:
# 1:4   -> Substraten
const Sind = 1:4
# 5:8   -> Enzymen
const Eind = 5:8
# N 
const Nind = Eind[end]+1
const Dind = Nind+1
const Lind = Dind+1
const lind = Lind+1
# P (Vrije fagen)
const Pfind = lind+1
const Paind = Pfind+1
const MOIind = Paind+1
# Benz (Benzonase)
const Benzind = MOIind+1

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
    fbaModel::Any
    # Resultaten voor Naïeve cellen (N)
    mu_N::Float64
    q_N::Vector{Float64}

    # Resultaten voor Lysogene cellen (l)
    mu_l::Float64
    q_l::Vector{Float64}
    q_benz_l::Float64  # Alleen de lysogenen produceren benzonase

    # --- NIEUWE ADSORPTIE PARAMETERS (LADDER) ---
    k_attach::Float64        # Binding rate (L / gDW / h) - voorheen alfa_ads
    k_dettach::Float64       # Loskoppelingsrate (1/h)
    k_inject::Float64       # Injectierate (1/h) - de stap van S_i naar I
    K_mal::Float64       # Affiniteit voor LamB (verzadiging van de ladder)


    infection_time::Float64

    # benzonase parameters
    benz_id::String    # Het ID van de reactie in de FBA
    k_tox::Float64     # Toxiciteitscoëfficiënt
    beta_benz::Float64 # Afbraaksnelheid enzym
    q_benz::Float64    # Opslag voor de berekende flux

    mu_max::Vector{Float64}   # [1.33, 1.26, 1.10, 0.29]
    e_max::Vector{Float64}    # Wordt berekend bij start
end


# 2. Receptor-factor (Storms: Michaelis-Menten verzadiging van LamB)
#function f_receptor(u, p::Parameters)
    # We halen het relatieve niveau van maltose-enzymen op
    #e_mal = u[p.ind_e[2]] 
    #return e_mal / (p.K_mal + e_mal + 1e-15)
#end


## --- POPULATIE DYNAMICA ---

function getPhi(u, p::Parameters)
    return u[MOIind]
end

function getProbLys(u, p::Parameters)
    phi = getPhi(u, p)
    return 1 - exp(-phi) - (phi * exp(-phi))
end

function getTotalBiomass(u, p::Parameters)
    return u[Nind] + u[Dind] + u[Lind] + u[lind]
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


function getPhages(u, p::Parameters)
    return u[Pfind]
end

