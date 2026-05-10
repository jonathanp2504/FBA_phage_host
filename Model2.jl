module Model2
import OrdinaryDiffEqCore
if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end
using JuMP
using OrdinaryDiffEq
using DelayDiffEq
using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels

# ============================================================
#  State vector indices (23 elementen)
# ============================================================
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
    duration::Float64
    startingBiomass::Float64
    alpha_syn::Float64
    beta_deg::Float64
    K_s::Vector{Float64}
    V_max::Vector{Float64}
    p_pref::Vector{Float64}
    tau::Float64
    b::Float64
    E_coli_cellDW::Float64
    MW::Vector{Float64}
    h_release::Float64
    biomass_id::String
    ex_ids::Vector{String}
    all_exchanges::Vector{String}
    essentials::Vector{String}
    fbaModelNaive::FbaCache
    fbaModelLysogen::FbaCache
    mu_N::Float64
    q_N::Vector{Float64}
    mu_l::Float64
    q_l::Vector{Float64}
    q_benz_l::Float64
    k_attach::Float64
    k_dettach::Float64
    k_inject::Float64
    K_mal::Float64
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

struct FbaSolution
    fluxes::Dict{String, Float64}
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

function getReceptorFactor(u, p::Parameters)
    e_mal = u[Eind[2]]; e_mal_max = p.e_max[2]
    return (e_mal + 1e-4) / e_mal_max
end

function getMonod(substrates::Vector{Float64}, p::Parameters)::Vector{Float64}
    return [max(0.0, substrates[i] / (substrates[i] + p.K_s[i])) for i in eachindex(substrates)]
end

function getRate(monod::Vector{Float64}, p::Parameters)::Vector{Float64}
    return p.V_max .* monod
end

function getU_cyt(rate::Vector{Float64}, p::Parameters)::Vector{Float64}
    weighted_R = p.p_pref .* rate
    return weighted_R ./ (sum(weighted_R) + 1e-10)
end

function getV_cyt(rate::Vector{Float64})::Vector{Float64}
    return rate ./ (maximum(rate) + 1e-10)
end

function getFluxes(sol, ids::Vector{String})::Vector{Float64}
    return [haskey(sol.fluxes, id) ? sol.fluxes[id] : 0.0 for id in ids]
end

function loadFBAmodel(path)
    model = convert(AbstractFBCModels.CanonicalModel.Model, load_model(path))
    model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0
    return model
end

function buildFbaCache(model, exchange_ids::Vector{String}, biomass_id::String; benz_id::Union{Nothing,String}=nothing)
    optimizer = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(optimizer)
    reaction_vars = Dict{String, JuMP.VariableRef}()
    for (rid, rxn) in model.reactions
        reaction_vars[rid] = JuMP.@variable(optimizer, lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound, base_name=rid)
    end
    balances = Dict(mid => JuMP.AffExpr() for mid in keys(model.metabolites))
    for (rid, rxn) in model.reactions
        for (mid, coef) in rxn.stoichiometry
            JuMP.add_to_expression!(balances[mid], coef, reaction_vars[rid])
        end
    end
    for bal in values(balances); JuMP.@constraint(optimizer, bal == 0.0); end
    JuMP.@objective(optimizer, Max, reaction_vars[biomass_id])
    tracked_ids = copy(exchange_ids)
    push!(tracked_ids, biomass_id)
    if !isnothing(benz_id); push!(tracked_ids, benz_id); end
    unique!(tracked_ids)
    return FbaCache(optimizer, reaction_vars, copy(exchange_ids), tracked_ids, biomass_id, benz_id)
end

function setLowerBound!(cache::FbaCache, rid::String, val::Float64)
    JuMP.set_lower_bound(cache.reaction_vars[rid], val)
end

function solveFba(cache::FbaCache, vmax::Vector{Float64}, p::Parameters)
    for i in eachindex(vmax); setLowerBound!(cache, cache.exchange_ids[i], -vmax[i]); end
    JuMP.optimize!(cache.optimizer)
    !JuMP.is_solved_and_feasible(cache.optimizer; dual=false) && return nothing
    fluxes = Dict{String,Float64}(rid => JuMP.value(cache.reaction_vars[rid]) for rid in cache.tracked_ids)
    return FbaSolution(fluxes)
end

function getVmax(u, p::Parameters)
    substrates = u[Sind]; enzymes = u[Eind]
    f = getMonod(substrates, p); R = getRate(f, p); v_cyt = getV_cyt(R)
    vmax = zeros(length(substrates))
    for i in eachindex(substrates)
        vmax[i] = R[i] * (enzymes[i] / p.e_max[i]) * v_cyt[i]
    end
    return vmax
end

function fbaUpdate!(u::Vector{Float64}, p::Parameters)
    vmax = getVmax(u, p)
    sol = solveFba(p.fbaModelNaive, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_N = sol.fluxes[p.biomass_id]; p.q_N = getFluxes(sol, p.ex_ids)
    else
        p.mu_N = 0.0; p.q_N .= 0.0
    end
    sol = solveFba(p.fbaModelLysogen, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_l = sol.fluxes[p.biomass_id]; p.q_l = getFluxes(sol, p.ex_ids)
    else
        p.mu_l = 0.0; p.q_l .= 0.0
    end
end

function enforcePositiveDomain!(u)
    for i in eachindex(u); u[i] = max(0.0, u[i]); end
end

function run(p::Parameters)
    u0 = zeros(23)
    u0[Sind] = [4.44, 2.337, 5.42, 0.0]
    u0[Eind] = [0.95, 0.01, 0.01, 0.01]
    u0[Nind] = p.startingBiomass
    tspan = (0.0, p.duration)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)               = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator) = any(x -> x < 0.0, u)
    domainAffect!(integrator)         = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(simulate_dFBA!, u0, (p, t) -> u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=[p.infection_time; fbaUpdateTimepoints],
        callback=CallbackSet(domainCallBack, infectionCallBack, fbaCallBack))
end

function simulate_dFBA!(du, u, h, p::Parameters, t)::Nothing
    updatePhageHostRates!(du, u, h, p, t)
    updateSubstrateUptakeRates!(du, u, h, p, t)
    return nothing
end

function updatePhageHostRates!(du, u, h, p::Parameters, t)::Nothing
    uDecision = h(p, t - 20/60)
    uLysis    = h(p, t - 60/60)
    X_tot     = getTotalBiomass(u, p)
    f_receptor = getReceptorFactor(u, p)
    k_attach_eff = p.k_attach * f_receptor
    k_inject_eff = p.k_inject * f_receptor
    mu_eff_l = p.mu_l * (1.0 - p.f_prod)

    du[Nind]  = p.mu_N * u[Nind] - k_inject_eff * u[Paind] * u[Nind] / X_tot
    du[Dind]  = k_inject_eff * u[Paind] * u[Nind] / X_tot
    du[Dind] -= k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[Lind]  = (1 - getProbLys(u, p)) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[Lind] -= (1 - getProbLys(uDecision, p)) * k_inject_eff * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    du[lind]  = getProbLys(u, p) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[lind] += mu_eff_l * u[lind] - p.k_tox * u[lind]
    du[MOIind]  = k_inject_eff * u[Paind] / X_tot
    du[MOIind] -= k_inject_eff * uDecision[Paind] / getTotalBiomass(uDecision, p)
    du[Pfind]  = p.b * (1 - getProbLys(uDecision, p)) * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    du[Pfind] -= k_attach_eff * X_tot * u[Pfind]
    du[Pfind] += p.k_dettach * u[Paind]
    du[Paind]  = k_attach_eff * X_tot * u[Pfind] - p.k_dettach * u[Paind] - k_inject_eff * u[Paind]
    du[Sind[1]] = p.h_release * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    groeiverlies = p.mu_l - mu_eff_l
    du[Benzind]  = groeiverlies * u[lind] * p.Y_benz
    return nothing
end

function updateSubstrateUptakeRates!(du, u, h, p::Parameters, t)::Nothing
    substrates = u[Sind]; enzymes = u[Eind]
    f = getMonod(substrates, p); R = getRate(f, p); u_cyt = getU_cyt(R, p)
    X_tot = getTotalBiomass(u, p)
    for i in eachindex(substrates)
        sub_idx = Sind[i]; enz_idx = Eind[i]
        flux_totaal = p.q_N[i] * X_tot * p.E_coli_cellDW + p.q_l[i] * u[lind] * p.E_coli_cellDW
        du[sub_idx] = min(0.0, flux_totaal)
        mu_avg = X_tot > 1.0 ? (p.mu_N * u[Nind] + p.mu_l * u[lind]) / X_tot : p.mu_N
        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - (p.beta_deg + mu_avg) * enzymes[i] + 0.001
    end
    return nothing
end

end # module Model2
