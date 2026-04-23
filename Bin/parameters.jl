using JuMP

# ============================================================
#  State vector indices
#  [Glc, Mal, Glyc, Ac,       -> Sind  1:4
#   e_glc, e_mal, e_gly, e_ac -> Eind  5:8
#   S, I, L,                  -> Sind  9,10,11
#   P,                        -> Pfind 12
#   Benz]                     -> Benzind 13
# ============================================================
const Sind    = 1:4
const Eind    = 5:8
const Sind_S  = Eind[end]+1   # 9  Susceptible (naïef)
const Sind_I  = Sind_S+1      # 10 Infected (lytisch)
const Sind_L  = Sind_I+1      # 11 Lysogeen
const Pfind   = Sind_L+1      # 12 Vrije fagen
const Benzind = Pfind+1       # 13 Benzonase

mutable struct FbaCache
    optimizer::JuMP.Model
    reaction_vars::Dict{String, JuMP.VariableRef}
    exchange_ids::Vector{String}
    tracked_ids::Vector{String}
    biomass_id::String
    benz_id::Union{Nothing, String}
end

mutable struct Parameters
    # Simulatie instellingen
    duration::Float64
    startingBiomass::Float64

    # Cybernetica & Kinetiek
    alpha_syn::Float64
    beta_deg::Float64
    K_s::Vector{Float64}
    V_max::Vector{Float64}
    p_pref::Vector{Float64}

    # Phage-Host Parameters (one-step adsorption)
    tau::Float64          # Latente periode (uren)
    b::Float64            # Burst size
    alfa_ads::Float64     # Adsorptieconstante (L/cel/h)
    E_coli_cellDW::Float64
    MW::Vector{Float64}
    h_release::Float64

    # Model IDs
    biomass_id::String
    ex_ids::Vector{String}
    all_exchanges::Vector{String}
    essentials::Vector{String}
    fbaModel::FbaCache

    # FBA resultaten
    mu::Float64
    q::Vector{Float64}

    infection_time::Float64
    infection_dose::Float64

    # Benzonase parameters (productie buiten FBA)
    benz_id::String
    k_tox::Float64
    beta_benz::Float64
    q_benz::Float64

    mu_max::Vector{Float64}
    e_max::Vector{Float64}

    # Groei-verlies methode voor Benzonase
    f_prod::Float64   # Fractie van groei opgeofferd aan Benzonase
    Y_benz::Float64   # mmol Benzonase per gDW groei-verlies
end

# ============================================================
#  Helper functies
# ============================================================

function getPhi(u, p::Parameters)
    alfa_eff = getAlfa_eff(u, p)
    return (alfa_eff * u[Pfind]) * 0.1
end

function getProbLys(phi::Float64)
    return 1 - exp(-phi) - (phi * exp(-phi))
end

function getTotalBiomass(u, p::Parameters)
    return u[Sind_S] + u[Sind_I] + u[Sind_L]
end

function getPhages(u, p::Parameters)
    return u[Pfind]
end

function getTotalAdsorptionFlux(u, p::Parameters)
    X_tot    = getTotalBiomass(u, p)
    alfa_eff = getAlfa_eff(u, p)
    return alfa_eff * X_tot * u[Pfind]
end

function getNewInfectionFlux(u, p::Parameters)
    alfa_eff = getAlfa_eff(u, p)
    return alfa_eff * u[Sind_S] * u[Pfind]
end

function getAlfa_eff(u, p::Parameters)
    e_mal     = u[Eind[2]]
    base_leak = 0.001
    return p.alfa_ads * (base_leak + (1 - base_leak) * e_mal)
end

function getMonod(substrates::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return [max(0.0, substrates[i] / (substrates[i] + parameters.K_s[i])) for i in eachindex(substrates)]
end

function getRate(monod::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return parameters.V_max .* monod
end

function getU_cyt(rate::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    weighted_R = parameters.p_pref .* rate
    denom = sum(weighted_R) + 1e-10
    return weighted_R ./ denom
end

function getV_cyt(rate::Vector{Float64})::Vector{Float64}
    denom = maximum(rate) + 1e-10
    return rate ./ denom
end

function getFluxes(sol, ids::Vector{String})::Vector{Float64}
    flux_vector = zeros(length(ids))
    for i in eachindex(ids)
        id = ids[i]
        flux_vector[i] = haskey(sol.fluxes, id) ? sol.fluxes[id] : 0.0
    end
    return flux_vector
end