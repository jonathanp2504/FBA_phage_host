module DEVerbosityPatch
    using OrdinaryDiffEqCore
    if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
        const DEVerbosity = false
        export DEVerbosity
    end
end
using .DEVerbosityPatch

include("./parameters.jl")
include("./FBA.jl")
using UnPack
using OrdinaryDiffEq
import OrdinaryDiffEqCore
using DelayDiffEq

# ============================================================
#  run() — hoofdfunctie, zelfde structuur als codeRefactorDries
# ============================================================
function run(p::Parameters)
    u0          = zeros(13)
    u0[Sind]    = [4.44, 2.337, 5.42, 0.0]
    u0[Eind]    = [0.95, 0.01, 0.01, 0.01]
    u0[Sind_S]  = p.startingBiomass

    tspan = (0.0, p.duration)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints              = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)           = fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator) = any(x -> x < 0.0, u)
    domainAffect!(integrator)         = enforcePositiveDomain!(integrator.u)
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

# ============================================================
#  ODE rechterkant
# ============================================================
function simulate_dFBA!(du, u, h, p::Parameters, t)::Nothing
    S_subs = u[Sind]
    e_enz  = u[Eind]
    S_cell = u[Sind_S]
    I_cell = u[Sind_I]
    L_cell = u[Sind_L]
    C_benz = u[Benzind]

    f     = getMonod(S_subs, p)
    R     = getRate(f, p)
    u_cyt = getU_cyt(R, p)

    phi      = getPhi(u, p)
    prob_lys = getProbLys(phi)

    X_tot                = getTotalBiomass(u, p)
    totaal_adsorptie     = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie      = getNewInfectionFlux(u, p)

    u_p        = h(p, t - p.tau)
    lysis_term = t > p.tau ? getNewInfectionFlux(u_p, p) : 0.0

    mu_safe    = max(0.0, p.mu)

    # --- Benzonase groei-verlies methode ---
    mu_eff_l       = mu_safe * (1.0 - p.f_prod)
    groeiverlies   = mu_safe - mu_eff_l
    productie_benz = groeiverlies * L_cell * p.Y_benz

    # --- Substraten & enzymen ---
    for i in eachindex(Sind)
        sub_idx = Sind[i]
        enz_idx = Eind[i]

        if S_subs[i] < 1e-8 && p.q[i] >= 0
            du[sub_idx] = 0.0
        else
            du[sub_idx] = p.q[i] * p.E_coli_cellDW * X_tot
        end

        if i == 1
            du[sub_idx] += p.h_release * lysis_term
        end

        if X_tot > 1.0
            mu_avg = (mu_safe * S_cell + mu_eff_l * L_cell) / X_tot
        else
            mu_avg = mu_safe
        end

        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] - (p.beta_deg + mu_avg) * e_enz[i] + 0.001
    end

    # --- Populatiedynamica (one-step adsorption) ---
    du[Sind_S] = mu_safe  * S_cell - nieuwe_infectie
    du[Sind_I] = (1 - prob_lys) * nieuwe_infectie - lysis_term
    du[Sind_L] = mu_eff_l * L_cell + prob_lys * nieuwe_infectie

    du[Pfind]   = p.b * lysis_term - totaal_adsorptie

    # --- Benzonase ---
    du[Benzind] = productie_benz

    return nothing
end

function enforcePositiveDomain!(u)
    for i in eachindex(u)
        u[i] = max(0.0, u[i])
    end
end