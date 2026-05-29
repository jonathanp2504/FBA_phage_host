module Model1
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
#  State vector indices (13 elementen)
# ============================================================
const Sind    = 1:4
const Eind    = 5:8
const Sind_S  = Eind[end]+1   # 9
const Sind_I  = Sind_S+1      # 10
const Sind_L  = Sind_I+1      # 11
const Pfind   = Sind_L+1      # 12
const Benzind = Pfind+1       # 13

const CELL_THRESHOLD  = 100.0
const PHAGE_THRESHOLD = 100.0

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

struct FbaSolution
    fluxes::Dict{String, Float64}
end

function getPhi(u, p::Parameters)
    S = u[Sind_S]; P = u[Pfind]
    return S > 1.0 ? min(P / S, 10.0) : 0.0
end

function getProbLys(phi::Float64)
    phi = max(0.0, phi)
    p = 1.0 - exp(-phi) - (phi * exp(-phi))
    return p < 0.01 ? 0.0 : p
end

function getTotalBiomass(u, p::Parameters)
    return u[Sind_S] + u[Sind_I] + u[Sind_L]
end

function getTotalAdsorptionFlux(u, p::Parameters)
    return getAlfa_eff(u, p) * getTotalBiomass(u, p) * u[Pfind]
end

function getNewInfectionFlux(u, p::Parameters)
    return getAlfa_eff(u, p) * u[Sind_S] * u[Pfind]
end

function getAlfa_eff(u, p::Parameters)
    e_mal     = u[Eind[2]]
    e_mal_max = p.e_max[2]
    f_rec     = e_mal / e_mal_max  # genormaliseerd zoals Model 2
    return p.alfa_ads * f_rec
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
    for bal in values(balances)
        JuMP.@constraint(optimizer, bal == 0.0)
    end
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
    for i in eachindex(vmax)
        setLowerBound!(cache, cache.exchange_ids[i], -vmax[i])
    end
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
    sol = solveFba(p.fbaModel, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu = sol.fluxes[p.biomass_id]
        p.q  = getFluxes(sol, p.ex_ids)
    else
        p.mu = 0.0; p.q .= 0.0
    end
    p.q_benz = 0.0
end

function enforcePositiveDomain!(u)
    for i in eachindex(u); u[i] = max(0.0, u[i]); end
end

function enforceThreshold!(u)
    if u[Sind_S] < CELL_THRESHOLD   u[Sind_S] = 0.0 end
    if u[Sind_I] < CELL_THRESHOLD   u[Sind_I] = 0.0 end
    if u[Sind_L] < CELL_THRESHOLD   u[Sind_L] = 0.0 end
    if u[Pfind]  < PHAGE_THRESHOLD  u[Pfind]  = 0.0 end
end

function run(p::Parameters)
    u0 = zeros(13)
    u0[Sind] = [4.44, 2.337, 5.42, 0.0]
    u0[Eind] = [0.95, 0.01, 0.01, 0.01]
    u0[Sind_S] = p.startingBiomass
    tspan = (0.0, p.duration)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
    function fbaAffect!(integrator)
        fbaUpdate!(integrator.u, p)
        enforceThreshold!(integrator.u)
    end
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
    S_subs = u[Sind]; e_enz = u[Eind]
    S_cell = u[Sind_S]; I_cell = u[Sind_I]; L_cell = u[Sind_L]
    f = getMonod(S_subs, p); R = getRate(f, p); u_cyt = getU_cyt(R, p)

    u_tau = h(p, t - p.tau)
    S_tau = u_tau[Sind_S]; P_tau = u_tau[Pfind]
    lysis_term = (t > p.tau && S_tau > CELL_THRESHOLD && P_tau > PHAGE_THRESHOLD) ?
        getNewInfectionFlux(u_tau, p) : 0.0
    phi_beslissing = (t > p.tau && S_tau > CELL_THRESHOLD) ? getPhi(u_tau, p) : 0.0
    prob_lys = getProbLys(phi_beslissing)

    X_tot            = getTotalBiomass(u, p)
    totaal_adsorptie = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie  = getNewInfectionFlux(u, p)
    mu_safe  = max(0.0, p.mu)
    mu_eff_l = mu_safe * (1.0 - p.f_prod)
    productie_benz = (mu_safe - mu_eff_l) * L_cell * p.Y_benz

    for i in eachindex(Sind)
        sub_idx = Sind[i]; enz_idx = Eind[i]
        du[sub_idx] = S_subs[i] < 1e-8 && p.q[i] >= 0 ? 0.0 : p.q[i] * p.E_coli_cellDW * X_tot
        if i == 1; du[sub_idx] += p.h_release * lysis_term; end
        mu_avg = X_tot > 1.0 ? (mu_safe * S_cell + mu_eff_l * L_cell) / X_tot : mu_safe
        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - (p.beta_deg + mu_avg) * e_enz[i] + 0.001
    end

    S_actief = S_cell > CELL_THRESHOLD ? S_cell : 0.0
    L_actief = L_cell > CELL_THRESHOLD ? L_cell : 0.0
    du[Sind_S]  = S_actief > 0.0 ? mu_safe * S_actief - nieuwe_infectie : 0.0
    du[Sind_I]  = (1.0 - prob_lys) * nieuwe_infectie - lysis_term
    du[Sind_L]  = mu_eff_l * L_actief + prob_lys * nieuwe_infectie
    du[Pfind]   = p.b * lysis_term - totaal_adsorptie
    du[Benzind] = productie_benz
    return nothing
end

end # module Model1
