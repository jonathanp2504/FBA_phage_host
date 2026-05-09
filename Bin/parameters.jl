using JuMP

# ============================================================
#  State vector indices (13 elementen)
#  [Glc, Mal, Glyc, Ac,       -> Sind    1:4
#   e_glc, e_mal, e_gly, e_ac -> Eind    5:8
#   S,                        -> Sind_S  9  Vatbaar
#   I,                        -> Sind_I  10 Lytisch
#   L,                        -> Sind_L  11 Lysogeen
#   P,                        -> Pfind   12 Vrije fagen
#   Benz]                     -> Benzind 13 Benzonase
#
#  MOI zit NIET in de toestandsvector.
#  Bereken post-hoc in de main als:
#    MOI = [sol.u[i][Pfind] / max(sol.u[i][Sind_S], 1.0) for i in eachindex(sol.u)]
# ============================================================
const Sind    = 1:4
const Eind    = 5:8
const Sind_S  = Eind[end]+1   # 9
const Sind_I  = Sind_S+1      # 10
const Sind_L  = Sind_I+1      # 11
const Pfind   = Sind_L+1      # 12
const Benzind = Pfind+1       # 13

mutable struct FbaCache
    optimizer::JuMP.Model
    reaction_vars::Dict{String, JuMP.VariableRef}
    exchange_ids::Vector{String}
    tracked_ids::Vector{String}
    biomass_id::String
    benz_id::Union{Nothing, String}
end

mutable struct Parameters
    duration::Float64
    startingBiomass::Float64
    alpha_syn::Float64
    beta_deg::Float64
    K_s::Vector{Float64}
    V_max::Vector{Float64}
    p_pref::Vector{Float64}
    tau::Float64
    b::Float64
    alfa_ads::Float64
    E_coli_cellDW::Float64
    MW::Vector{Float64}
    h_release::Float64
    biomass_id::String
    ex_ids::Vector{String}
    all_exchanges::Vector{String}
    essentials::Vector{String}
    fbaModel::FbaCache
    mu::Float64
    q::Vector{Float64}
    infection_time::Float64
    infection_dose::Float64
    benz_id::String
    k_tox::Float64
    beta_benz::Float64
    q_benz::Float64
    mu_max::Vector{Float64}
    e_max::Vector{Float64}
    f_prod::Float64
    Y_benz::Float64
end

# ============================================================
#  Helper functies
# ============================================================

# MOI algebraïsch: P/S op dit moment.
# Cap op 10: getProbLys(phi>5) ≈ 1 toch al, hogere waarden
# veranderen de biologie niet maar veroorzaken numerieke problemen.
function getPhi(u, p::Parameters)
    S = u[Sind_S]
    P = u[Pfind]
    return S > 1.0 ? min(P / S, 10.0) : 0.0
end

# Kans op lysogenie als functie van MOI (phi).
# Gebaseerd op Poisson: P(>=2 infecties) = 1 - P(0) - P(1)
# phi << 1 -> prob_lys ≈ 0  -> lysis domineert
# phi >> 1 -> prob_lys ≈ 1  -> lysogenie domineert
function getProbLys(phi::Float64)
    phi = max(0.0, phi)
    p = 1.0 - exp(-phi) - (phi * exp(-phi))
    # Onder 1% kans op lysogenie -> pure lysis
    # Dit voorkomt numerieke accumulatie van microscopische L-waarden
    return p < 0.01 ? 0.0 : p
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
    return [max(0.0, substrates[i] / (substrates[i] + parameters.K_s[i]))
            for i in eachindex(substrates)]
end

function getRate(monod::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return parameters.V_max .* monod
end

function getU_cyt(rate::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    weighted_R = parameters.p_pref .* rate
    return weighted_R ./ (sum(weighted_R) + 1e-10)
end

function getV_cyt(rate::Vector{Float64})::Vector{Float64}
    return rate ./ (maximum(rate) + 1e-10)
end

function getFluxes(sol, ids::Vector{String})::Vector{Float64}
    return [haskey(sol.fluxes, id) ? sol.fluxes[id] : 0.0 for id in ids]
end