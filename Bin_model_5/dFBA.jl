import OrdinaryDiffEqCore
if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end
using OrdinaryDiffEq
using DelayDiffEq
using UnPack
include("./parameters.jl")
include("./FBA.jl")

function run(p::Parameters)
    u0        = zeros(16)
    u0[Sind]  = [0.0, 2.337, 0.0, 0.0]
    u0[Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Nind]  = p.startingBiomass
    tspan     = (0.0, p.duration)

    infectionCondition(u, t, integrator)  = t == p.infection_time
    infectionAffect!(integrator)          = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints                   = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator)  = t in fbaUpdateTimepoints
    fbaAffect!(integrator)                = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator)     = any(x -> x < 0.0, u)
    domainAffect!(integrator)             = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem  = DDEProblem(simulate_dFBA!, u0, (p, t) -> u0, tspan, p)
    solution = solve(problem, MethodOfSteps(Tsit5()),
        verbose  = false,
        reltol   = 1e-4,
        abstol   = 1e-6,
        tstops   = [p.infection_time; fbaUpdateTimepoints],
        callback = CallbackSet(domainCallBack, infectionCallBack, fbaCallBack))
    return solution
end

function simulate_dFBA!(du, u, h, p::Parameters, t)::Nothing
    updatePhageHostRates!(du, u, h, p, t)
    updateSubstrateUptakeRates!(du, u, h, p, t)
    return nothing
end

function updatePhageHostRates!(du, u, h, p::Parameters, t)::Nothing
    uDecision    = h(p, t - 28/60)
    uLysis       = h(p, t - 79/60)
    X_tot        = getTotalBiomass(u, p)
    f_receptor   = getReceptorFactor(u, p)
    k_attach_eff = p.k_attach * f_receptor
    k_inject_eff = p.k_inject * f_receptor

    # MODEL 5: burst size afgeleid uit lytische FBA
    # q_phage_L [mmol faagpartikels/gDW/h]
    # × E_coli_cellDW [gDW/cel]
    # × tau [h] (latente periode)
    # × (AVOGADRO/1000) [partikels/mmol]
    # = faagpartikels per cel per latente periode
    tau_eclipse = 28.0 / 60.0   # 28 minuten in uren
    tau_rise    = max(0.0, p.tau - tau_eclipse)
    b_fba = max(1.0, p.q_phage_L * p.E_coli_cellDW * tau_rise * (AVOGADRO * 1e-3))
 

    du[Nind]  = p.mu_N * u[Nind]
    du[Nind] -= k_inject_eff * u[Paind] * u[Nind] / X_tot

    du[Dind]  = k_inject_eff * u[Paind] * u[Nind] / X_tot
    du[Dind] -= k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)

    du[Lind]  = (1 - getProbLys(u, p)) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[Lind] -= (1 - getProbLys(uDecision, p)) * k_inject_eff * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)

    du[lind]  = getProbLys(u, p) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[lind] += p.mu_l * u[lind]
    du[lind] -= p.k_tox * u[lind]

    du[MOIind]  = k_inject_eff * u[Paind] / X_tot
    du[MOIind] -= k_inject_eff * uDecision[Paind] / getTotalBiomass(uDecision, p)

    # b_fba uit de lytische FBA
    du[Pfind]  = b_fba * (1 - getProbLys(uDecision, p)) * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    du[Pfind] -= k_attach_eff * X_tot * u[Pfind]
    du[Pfind] += p.k_dettach * u[Paind]

    du[Paind]  = k_attach_eff * X_tot * u[Pfind]
    du[Paind] -= p.k_dettach * u[Paind]
    du[Paind] -= k_inject_eff * u[Paind]

    du[Sind[1]] += p.h_release * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    return nothing
end

function updateSubstrateUptakeRates!(du, u, h, p::Parameters, t)::Nothing
    substrates = u[Sind]
    enzymes    = u[Eind]
    f          = getMonod(substrates, p)
    R          = getRate(f, p)
    u_cyt      = getU_cyt(R, p)
    X_tot      = getTotalBiomass(u, p)

    for i in eachindex(substrates)
        sub_idx     = Sind[i]
        enz_idx     = Eind[i]
        flux_N      = p.q_N[i] * X_tot * p.E_coli_cellDW
        flux_l      = p.q_l[i] * u[lind] * p.E_coli_cellDW
        du[sub_idx] = min(0.0, flux_N + flux_l)
        mu_avg      = X_tot > 1.0 ? (p.mu_N * u[Nind] + p.mu_l * u[lind]) / X_tot : p.mu_N
        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - (p.beta_deg + mu_avg) * enzymes[i] + 0.001
    end

    du[Benzind] = p.q_benz_l * p.E_coli_cellDW * u[lind]
    return nothing
end

function getMonod(substrates::Vector{Float64}, parameters::Parameters)::Vector{Float64}
    return [max(0.0, substrates[i] / (substrates[i] + parameters.K_s[i])) for i in eachindex(substrates)]
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

function getReceptorFactor(u, p::Parameters)
    return (u[Eind[2]] + 1e-4) / p.e_max[2]
end

function enforcePositiveDomain!(u)
    for i in eachindex(u)
        u[i] = max(0.0, u[i])
    end
end
