using JuMP

# indices in state vector u:
const Sind    = 1:4
const Eind    = 5:8
const Nind    = Eind[end]+1
const Dind    = Nind+1
const Lind    = Dind+1
const lind    = Lind+1
const Pfind   = lind+1
const Paind   = Pfind+1
const MOIind  = Paind+1
const Benzind = MOIind+1

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

    # Phage-Host Parameters
    tau::Float64
    b::Float64
    E_coli_cellDW::Float64
    MW::Vector{Float64}
    h_release::Float64

    # Model IDs
    biomass_id::String
    ex_ids::Vector{String}
    all_exchanges::Vector{String}
    essentials::Vector{String}
    fbaModelNaive::FbaCache
    fbaModelLysogen::FbaCache

    # FBA resultaten naïeve cellen
    mu_N::Float64
    q_N::Vector{Float64}

    # FBA resultaten lysogene cellen
    mu_l::Float64
    q_l::Vector{Float64}
    q_benz_l::Float64

    # Adsorptie parameters
    k_attach::Float64
    k_dettach::Float64
    k_inject::Float64
    K_mal::Float64

    infection_time::Float64
    infection_dose::Float64

    # Benzonase parameters
    benz_id::String
    k_tox::Float64
    beta_benz::Float64
    q_benz::Float64

    mu_max::Vector{Float64}
    e_max::Vector{Float64}

    # Benzonase productie via groei-verlies methode (main branch aanpak)
    # Benzonase zit NIET in de FBA — productie wordt berekend in de ODE
    f_prod::Float64   # Fractie van groei opgeofferd aan Benzonase (bijv. 0.2)
    Y_benz::Float64   # Opbrengst: mmol Benzonase per gDW groei-verlies
end

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

function getPhages(u, p::Parameters)
    return u[Pfind]
end

function getMu_avg(p, u)
    total_cells = u[Nind] + u[lind]
    if total_cells > 1e-6
        return (p.mu_N * u[Nind] + p.mu_l * u[lind]) / total_cells
    else
        return p.mu_N
    end
end

benz_stoich = Dict(
    "M_ala__L_c" => -33.0,
    "M_arg__L_c" => -13.0,
    "M_asn__L_c" => -22.0,
    "M_asp__L_c" => -16.0,
    "M_cys__L_c" => -4.0,
    "M_gln__L_c" => -12.0,
    "M_glu__L_c" => -10.0,
    "M_gly_c"    => -21.0,
    "M_his__L_c" => -4.0,
    "M_ile__L_c" => -8.0,
    "M_leu__L_c" => -21.0,
    "M_lys__L_c" => -14.0,
    "M_met__L_c" => -3.0,
    "M_phe__L_c" => -8.0,
    "M_pro__L_c" => -10.0,
    "M_ser__L_c" => -19.0,
    "M_thr__L_c" => -15.0,
    "M_trp__L_c" => -5.0,
    "M_tyr__L_c" => -10.0,
    "M_val__L_c" => -10.0,
    "M_atp_c"    => -532.0,
    "M_adp_c"    =>  532.0,
    "M_gtp_c"    => -532.0,
    "M_gdp_c"    =>  532.0,
    "M_h2o_c"    => -1064.0,
    "M_pi_c"     =>  1064.0,
    "M_h_c"      =>  1064.0,
    "M_benzonase_c" => 1.0
)