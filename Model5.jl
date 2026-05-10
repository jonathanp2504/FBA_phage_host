module Model5
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
#  State vector: 16 elementen
# ============================================================
const Sind    = 1:4
const Eind    = 5:8
const Nind    = Eind[end]+1   # 9
const Dind    = Nind+1        # 10
const Lind    = Dind+1        # 11
const lind    = Lind+1        # 12
const Pfind   = lind+1        # 13
const Paind   = Pfind+1       # 14
const MOIind  = Paind+1       # 15
const Benzind = MOIind+1      # 16

const AVOGADRO = 6.022e23

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
    E_coli_cellDW::Float64
    MW::Vector{Float64}
    h_release::Float64
    biomass_id::String
    ex_ids::Vector{String}
    all_exchanges::Vector{String}
    essentials::Vector{String}
    fbaModelNaive::FbaCache
    fbaModelLysogen::FbaCache
    fbaModelLytic::FbaCache
    mu_N::Float64
    q_N::Vector{Float64}
    mu_l::Float64
    q_l::Vector{Float64}
    q_benz_l::Float64
    q_phage_L::Float64
    k_attach::Float64
    k_dettach::Float64
    k_inject::Float64
    K_mal::Float64
    infection_time::Float64
    infection_dose::Float64
    benz_id::String
    phage_id::String
    k_tox::Float64
    beta_benz::Float64
    q_benz::Float64
    mu_max::Vector{Float64}
    e_max::Vector{Float64}
    f_prod::Float64
end

struct FbaSolution
    fluxes::Dict{String, Float64}
end

function getPhi(u, p::Parameters); return u[MOIind]; end
function getProbLys(u, p::Parameters)
    phi = getPhi(u, p)
    return 1 - exp(-phi) - (phi * exp(-phi))
end
function getTotalBiomass(u, p::Parameters)
    return u[Nind] + u[Dind] + u[Lind] + u[lind]
end
function getReceptorFactor(u, p::Parameters)
    return (u[Eind[2]] + 1e-4) / p.e_max[2]
end
function getMonod(s::Vector{Float64}, p::Parameters)::Vector{Float64}
    return [max(0.0, s[i]/(s[i]+p.K_s[i])) for i in eachindex(s)]
end
function getRate(m::Vector{Float64}, p::Parameters)::Vector{Float64}
    return p.V_max .* m
end
function getU_cyt(r::Vector{Float64}, p::Parameters)::Vector{Float64}
    w = p.p_pref .* r; return w ./ (sum(w)+1e-10)
end
function getV_cyt(r::Vector{Float64})::Vector{Float64}
    return r ./ (maximum(r)+1e-10)
end
function getFluxes(sol, ids::Vector{String})::Vector{Float64}
    return [haskey(sol.fluxes, id) ? sol.fluxes[id] : 0.0 for id in ids]
end

function loadFBAmodel(path)
    model = convert(AbstractFBCModels.CanonicalModel.Model, load_model(path))
    model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0
    return model
end

function addBenzonase!(model, benz_stoich)
    model.metabolites["M_benzonase_c"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_benzonase_c"].name = "Benzonase"
    model.metabolites["M_benzonase_c"].compartment = "c"
    model.metabolites["M_benzonase_e"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_benzonase_e"].name = "Benzonase (extracellular)"
    model.metabolites["M_benzonase_e"].compartment = "e"
    model.reactions["R_BENZ_prod"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_BENZ_prod"].name = "Benzonase production"
    model.reactions["R_BENZ_prod"].lower_bound = 0.0
    model.reactions["R_BENZ_prod"].upper_bound = 1000.0
    model.reactions["R_BENZ_prod"].stoichiometry = benz_stoich
    model.reactions["R_BENZ_export"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_BENZ_export"].stoichiometry = Dict("M_benzonase_c"=>-1.0,"M_benzonase_e"=>1.0)
    model.reactions["R_BENZ_export"].lower_bound = 0.0
    model.reactions["R_BENZ_export"].upper_bound = 1000.0
    model.reactions["R_EX_benz_e"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_EX_benz_e"].stoichiometry = Dict("M_benzonase_e"=>-1.0)
    model.reactions["R_EX_benz_e"].lower_bound = 0.0
    model.reactions["R_EX_benz_e"].upper_bound = 1000.0
    return model
end

function addPhage!(model, phage_stoich)
    model.metabolites["M_phage_c"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_phage_c"].name = "Phage lambda virion"
    model.metabolites["M_phage_c"].compartment = "c"
    model.metabolites["M_phage_e"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_phage_e"].name = "Phage lambda (extracellular)"
    model.metabolites["M_phage_e"].compartment = "e"
    model.reactions["R_PHAGE_prod"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_PHAGE_prod"].name = "Phage lambda assembly"
    model.reactions["R_PHAGE_prod"].lower_bound = 0.0
    model.reactions["R_PHAGE_prod"].upper_bound = 1000.0
    model.reactions["R_PHAGE_prod"].stoichiometry = phage_stoich
    model.reactions["R_PHAGE_release"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_PHAGE_release"].stoichiometry = Dict("M_phage_c"=>-1.0,"M_phage_e"=>1.0)
    model.reactions["R_PHAGE_release"].lower_bound = 0.0
    model.reactions["R_PHAGE_release"].upper_bound = 1000.0
    model.reactions["R_EX_phage_e"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_EX_phage_e"].stoichiometry = Dict("M_phage_e"=>-1.0)
    model.reactions["R_EX_phage_e"].lower_bound = 0.0
    model.reactions["R_EX_phage_e"].upper_bound = 1000.0
    return model
end

const benz_stoich = Dict(
    "M_ala__L_c"=>-33.0,"M_arg__L_c"=>-13.0,"M_asn__L_c"=>-22.0,
    "M_asp__L_c"=>-16.0,"M_cys__L_c"=>-4.0, "M_gln__L_c"=>-12.0,
    "M_glu__L_c"=>-10.0,"M_gly_c"   =>-21.0,"M_his__L_c"=>-4.0,
    "M_ile__L_c"=>-8.0, "M_leu__L_c"=>-21.0,"M_lys__L_c"=>-14.0,
    "M_met__L_c"=>-3.0, "M_phe__L_c"=>-8.0, "M_pro__L_c"=>-10.0,
    "M_ser__L_c"=>-19.0,"M_thr__L_c"=>-15.0,"M_trp__L_c"=>-5.0,
    "M_tyr__L_c"=>-10.0,"M_val__L_c"=>-10.0,
    "M_atp_c"=>-532.0,"M_adp_c"=>532.0,"M_gtp_c"=>-532.0,"M_gdp_c"=>532.0,
    "M_h2o_c"=>-1064.0,"M_pi_c"=>1064.0,"M_h_c"=>1064.0,"M_benzonase_c"=>1.0
)

const phage_stoich = Dict(
    "M_ala__L_c"=>-50625.0,"M_arg__L_c"=>-8100.0,"M_asn__L_c"=>-6885.0,
    "M_asp__L_c"=>-6885.0, "M_gln__L_c"=>-36045.0,"M_glu__L_c"=>-8100.0,
    "M_gly_c"   =>-7290.0, "M_his__L_c"=>-405.0,  "M_ile__L_c"=>-10530.0,
    "M_leu__L_c"=>-20655.0,"M_lys__L_c"=>-12150.0,"M_met__L_c"=>-2835.0,
    "M_phe__L_c"=>-3240.0, "M_pro__L_c"=>-3240.0, "M_ser__L_c"=>-13365.0,
    "M_thr__L_c"=>-14985.0,"M_tyr__L_c"=>-5265.0, "M_val__L_c"=>-8100.0,
    "M_atp_c"=>-437400.0,"M_adp_c"=>437400.0,"M_gtp_c"=>-437400.0,"M_gdp_c"=>437400.0,
    "M_h2o_c"=>-874800.0,"M_pi_c"=>874800.0,"M_h_c"=>874800.0,
    "M_datp_c"=>-12149.0,"M_dttp_c"=>-12149.0,"M_dgtp_c"=>-12101.0,"M_dctp_c"=>-12101.0,
    "M_phage_c"=>1.0
)

function buildFbaCache(model, exchange_ids::Vector{String}, biomass_id::String;
                       benz_id::Union{Nothing,String}=nothing)
    optimizer = JuMP.Model(HiGHS.Optimizer); JuMP.set_silent(optimizer)
    reaction_vars = Dict{String,JuMP.VariableRef}()
    for (rid,rxn) in model.reactions
        reaction_vars[rid] = JuMP.@variable(optimizer,
            lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound, base_name=rid)
    end
    balances = Dict(mid=>JuMP.AffExpr() for mid in keys(model.metabolites))
    for (rid,rxn) in model.reactions
        for (mid,coef) in rxn.stoichiometry
            JuMP.add_to_expression!(balances[mid], coef, reaction_vars[rid])
        end
    end
    for bal in values(balances); JuMP.@constraint(optimizer, bal==0.0); end
    JuMP.@objective(optimizer, Max, reaction_vars[biomass_id])
    tracked_ids = copy(exchange_ids); push!(tracked_ids, biomass_id)
    if !isnothing(benz_id); push!(tracked_ids, benz_id); end
    unique!(tracked_ids)
    return FbaCache(optimizer, reaction_vars, copy(exchange_ids), tracked_ids, biomass_id, benz_id)
end

function buildLyticFbaCache(model, exchange_ids::Vector{String}, biomass_id::String, phage_id::String)
    optimizer = JuMP.Model(HiGHS.Optimizer); JuMP.set_silent(optimizer)
    reaction_vars = Dict{String,JuMP.VariableRef}()
    for (rid,rxn) in model.reactions
        reaction_vars[rid] = JuMP.@variable(optimizer,
            lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound, base_name=rid)
    end
    balances = Dict(mid=>JuMP.AffExpr() for mid in keys(model.metabolites))
    for (rid,rxn) in model.reactions
        for (mid,coef) in rxn.stoichiometry
            JuMP.add_to_expression!(balances[mid], coef, reaction_vars[rid])
        end
    end
    for bal in values(balances); JuMP.@constraint(optimizer, bal==0.0); end
    JuMP.@objective(optimizer, Max, reaction_vars[phage_id])
    tracked_ids = copy(exchange_ids)
    push!(tracked_ids, biomass_id); push!(tracked_ids, phage_id)
    unique!(tracked_ids)
    return FbaCache(optimizer, reaction_vars, copy(exchange_ids), tracked_ids, biomass_id, phage_id)
end

function setLowerBound!(cache::FbaCache, rid::String, val::Float64)
    JuMP.set_lower_bound(cache.reaction_vars[rid], val)
end

function solveFba(cache::FbaCache, vmax::Vector{Float64}, p::Parameters)
    for i in eachindex(vmax); setLowerBound!(cache, cache.exchange_ids[i], -vmax[i]); end
    JuMP.optimize!(cache.optimizer)
    !JuMP.is_solved_and_feasible(cache.optimizer; dual=false) && return nothing
    fluxes = Dict{String,Float64}(rid=>JuMP.value(cache.reaction_vars[rid]) for rid in cache.tracked_ids)
    return FbaSolution(fluxes)
end

function getVmax(u, p::Parameters)
    s=u[Sind]; e=u[Eind]; f=getMonod(s,p); R=getRate(f,p); v=getV_cyt(R)
    vmax=zeros(length(s))
    for i in eachindex(s); vmax[i]=R[i]*(e[i]/p.e_max[i])*v[i]; end
    return vmax
end

function fbaUpdate!(u::Vector{Float64}, p::Parameters)
    vmax = getVmax(u, p)
    totale_import = sum(vmax)
    sol = solveFba(p.fbaModelNaive, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_N = sol.fluxes[p.biomass_id]; p.q_N = getFluxes(sol, p.ex_ids)
    else
        p.mu_N = 0.0; p.q_N .= 0.0
    end
    setLowerBound!(p.fbaModelLysogen, p.benz_id, totale_import * p.f_prod)
    sol = solveFba(p.fbaModelLysogen, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_l = sol.fluxes[p.biomass_id]
        p.q_l  = getFluxes(sol, p.ex_ids)
        p.q_benz_l = sol.fluxes[p.benz_id]
    else
        p.mu_l = 0.0; p.q_l .= 0.0; p.q_benz_l = 0.0
    end
    sol = solveFba(p.fbaModelLytic, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.q_phage_L = sol.fluxes[p.phage_id]
    else
        p.q_phage_L = 0.0
    end
end

function enforcePositiveDomain!(u)
    for i in eachindex(u); u[i] = max(0.0, u[i]); end
end

function run(p::Parameters)
    u0 = zeros(16)
    u0[Sind] = [0.0, 2.337, 5.42, 0.0]
    u0[Eind] = [0.95, 0.01, 0.01, 0.01]
    u0[Nind] = p.startingBiomass
    tspan = (0.0, p.duration)

    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)
    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)
    domainCondition(u,t,integrator) = any(x->x<0.0, u)
    domainAffect!(integrator)       = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(simulate_dFBA!, u0, (p,t)->u0, tspan, p)
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
    uDecision    = h(p, t - 20/60)
    uLysis       = h(p, t - 60/60)
    X_tot        = getTotalBiomass(u, p)
    f_receptor   = getReceptorFactor(u, p)
    k_attach_eff = p.k_attach * f_receptor
    k_inject_eff = p.k_inject * f_receptor
    b_fba = max(1.0, p.q_phage_L * p.E_coli_cellDW * p.tau * (AVOGADRO * 1e-3))

    du[Nind]  = p.mu_N * u[Nind] - k_inject_eff * u[Paind] * u[Nind] / X_tot
    du[Dind]  = k_inject_eff * u[Paind] * u[Nind] / X_tot
    du[Dind] -= k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[Lind]  = (1-getProbLys(u,p)) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[Lind] -= (1-getProbLys(uDecision,p)) * k_inject_eff * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    du[lind]  = getProbLys(u,p) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[lind] += p.mu_l * u[lind] - p.k_tox * u[lind]
    du[MOIind]  = k_inject_eff * u[Paind] / X_tot
    du[MOIind] -= k_inject_eff * uDecision[Paind] / getTotalBiomass(uDecision, p)
    du[Pfind]  = b_fba * (1-getProbLys(uDecision,p)) * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    du[Pfind] -= k_attach_eff * X_tot * u[Pfind]
    du[Pfind] += p.k_dettach * u[Paind]
    du[Paind]  = k_attach_eff * X_tot * u[Pfind] - p.k_dettach * u[Paind] - k_inject_eff * u[Paind]
    du[Sind[1]] += p.h_release * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    return nothing
end

function updateSubstrateUptakeRates!(du, u, h, p::Parameters, t)::Nothing
    s=u[Sind]; e=u[Eind]; f=getMonod(s,p); R=getRate(f,p); u_cyt=getU_cyt(R,p)
    X_tot = getTotalBiomass(u, p)
    for i in eachindex(s)
        du[Sind[i]] = min(0.0, p.q_N[i]*X_tot*p.E_coli_cellDW + p.q_l[i]*u[lind]*p.E_coli_cellDW)
        mu_avg = X_tot>1.0 ? (p.mu_N*u[Nind]+p.mu_l*u[lind])/X_tot : p.mu_N
        du[Eind[i]] = p.alpha_syn*f[i]*u_cyt[i] - (p.beta_deg+mu_avg)*e[i] + 0.001
    end
    du[Benzind] = p.q_benz_l * p.E_coli_cellDW * u[lind]
    return nothing
end

end # module Model5