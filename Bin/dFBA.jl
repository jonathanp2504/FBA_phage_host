import OrdinaryDiffEqCore
if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end
using OrdinaryDiffEq
using DelayDiffEq

include("./parameters.jl")
include("./FBA.jl")
using UnPack
using OrdinaryDiffEq
import OrdinaryDiffEqCore
using DelayDiffEq

function run(p::Parameters)
    u0 = zeros(23)
    u0[Sind] = [4.44, 2.337, 5.42, 0.0]   # glucose AAN in deze branch
    u0[Eind] = [0.95, 0.01, 0.01, 0.01]
    u0[Nind] = p.startingBiomass

    tspan = (0.0, p.duration)

    infectionCondition(u, t, integrator)  = t == p.infection_time
    infectionAffect!(integrator)          = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator)  = t in fbaUpdateTimepoints
    fbaAffect!(integrator)                = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator)     = any(x -> x < 0.0, u)
    domainAffect!(integrator)             = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(simulate_dFBA!, u0, (p, t) -> u0, tspan, p)

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
    uDecision  = h(p, t - 20/60)
    uLysis     = h(p, t - 60/60)
    X_tot      = getTotalBiomass(u, p)
    f_receptor = getReceptorFactor(u, p)

    k_attach_eff = p.k_attach * f_receptor
    k_inject_eff = p.k_inject * f_receptor

    # Effectieve groei lysogeen: verminderd door Benzonase-productielast
    mu_eff_l = p.mu_l * (1.0 - p.f_prod)

    du[Nind]  = p.mu_N * u[Nind]
    du[Nind] -= k_inject_eff * u[Paind] * u[Nind] / X_tot

    du[Dind]  = k_inject_eff * u[Paind] * u[Nind] / X_tot
    du[Dind] -= k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)

    du[Lind]  = (1 - getProbLys(u, p)) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[Lind] -= (1 - getProbLys(uDecision, p)) * k_inject_eff * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)

    du[lind]  = getProbLys(u, p) * k_inject_eff * uDecision[Paind] * uDecision[Nind] / getTotalBiomass(uDecision, p)
    du[lind] += mu_eff_l * u[lind]
    du[lind] -= p.k_tox * u[lind]

    du[MOIind]  = k_inject_eff * u[Paind] / X_tot
    du[MOIind] -= k_inject_eff * uDecision[Paind] / getTotalBiomass(uDecision, p)

    du[Pfind]  = p.b * (1 - getProbLys(uDecision, p)) * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)
    du[Pfind] -= k_attach_eff * X_tot * u[Pfind]
    du[Pfind] += p.k_dettach * u[Paind]

    du[Paind]  = k_attach_eff * X_tot * u[Pfind]
    du[Paind] -= p.k_dettach * u[Paind]
    du[Paind] -= k_inject_eff * u[Paind]

    du[Sind[1]] += p.h_release * p.k_inject * uLysis[Paind] * uLysis[Nind] / getTotalBiomass(uLysis, p)

    # --- BENZONASE PRODUCTIE (buiten FBA, via groei-verlies methode) ---
    groeiverlies   = p.mu_l - mu_eff_l
    productie_benz = groeiverlies * u[lind] * p.Y_benz
    du[Benzind]    = productie_benz

    return nothing
end

function updateSubstrateUptakeRates!(du, u, h, p::Parameters, t)::Nothing
    substrates = u[Sind]
    enzymes    = u[Eind]
    f     = getMonod(substrates, p)
    R     = getRate(f, p)
    u_cyt = getU_cyt(R, p)
    X_tot = getTotalBiomass(u, p)

    for i in eachindex(substrates)
        sub_idx = Sind[i]
        enz_idx = Eind[i]

        flux_N      = p.q_N[i] * X_tot * p.E_coli_cellDW
        flux_l      = p.q_l[i] * u[lind] * p.E_coli_cellDW
        flux_totaal = flux_N + flux_l
        du[sub_idx] = min(0.0, flux_totaal)

        if X_tot > 1.0
            mu_avg = (p.mu_N * u[Nind] + p.mu_l * u[lind]) / X_tot
        else
            mu_avg = p.mu_N
        end

        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - (p.beta_deg + mu_avg) * enzymes[i] + 0.001
    end
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

function getReceptorFactor(u, p::Parameters)
    e_mal     = u[Eind[2]]
    e_mal_max = p.e_max[2]
    return (e_mal + 1e-4) / e_mal_max
end

function enforcePositiveDomain!(u)
    for i in eachindex(u)
        u[i] = max(0.0, u[i])
    end
end