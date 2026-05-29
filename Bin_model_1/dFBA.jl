import OrdinaryDiffEqCore
if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end
using OrdinaryDiffEq
using DelayDiffEq
using UnPack

include("./parameters.jl")
include("./FBA.jl")

const CELL_THRESHOLD  = 100.0
const PHAGE_THRESHOLD = 100.0

function run(p::Parameters)
    u0         = zeros(13)
    u0[Sind]   = [4.44, 2.337, 0.0, 0.0]
    u0[Eind]   = [0.95, 0.01, 0.01, 0.01]
    u0[Sind_S] = p.startingBiomass

    tspan = (0.0, p.duration)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints                  = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints

    function fbaAffect!(integrator)
        fbaUpdate!(integrator.u, p)
        enforceThreshold!(integrator.u)
    end

    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u, t, integrator)    = any(x -> x < 0.0, u)
    domainAffect!(integrator)            = enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem  = DDEProblem(simulate_dFBA!, u0, (p, t) -> u0, tspan, p)
    solution = solve(problem, MethodOfSteps(Tsit5()),
        verbose  = false,
        reltol   = 1e-4,
        abstol   = 1e-6,
        tstops = [p.infection_time;
          p.infection_time + (28/60);
          p.infection_time + p.tau;
          fbaUpdateTimepoints],
        callback = CallbackSet(domainCallBack, infectionCallBack, fbaCallBack))
    return solution
end

function simulate_dFBA!(du, u, h, p::Parameters, t)::Nothing
    S_subs = u[Sind]
    e_enz  = u[Eind]
    S_cell = u[Sind_S]
    I_cell = u[Sind_I]
    L_cell = u[Sind_L]

    f     = getMonod(S_subs, p)
    R     = getRate(f, p)
    u_cyt = getU_cyt(R, p)

    # Twee aparte vertraagde toestanden
    # En in simulate_dFBA! gebruik:
    tau_d = 28.0 / 60.0   # altijd vast op 28 minuten, biologische constante
    tau_l = p.tau          # varieert in sensitivity analyse  # 28 minuten in

    u_decision = h(p, t - tau_d)
    u_lysis    = h(p, t - tau_l)

    # Lysis term: gebaseerd op infectieflux van TAU_L geleden
    S_lysis = u_lysis[Sind_S]
    P_lysis = u_lysis[Pfind]
    # Door:
    lysis_term = (t > p.infection_time + tau_l && S_lysis > CELL_THRESHOLD && P_lysis > PHAGE_THRESHOLD) ?
        getNewInfectionFlux(u_lysis, p) : 0.0

    # Beslissing: gebaseerd op MOI van TAU_D geleden (moment van infectie)
    S_decision = u_decision[Sind_S]
    phi_beslissing = (t > p.infection_time + tau_d && S_decision > CELL_THRESHOLD) ?
        getPhi(u_decision, p) : 0.0
    prob_lys = getProbLys(phi_beslissing)

    X_tot            = getTotalBiomass(u, p)
    totaal_adsorptie = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie  = getNewInfectionFlux(u, p)

    mu_safe  = max(0.0, p.mu)
    mu_eff_l = mu_safe * (1.0 - p.f_prod)

    productie_benz = (mu_safe - mu_eff_l) * L_cell * p.Y_benz

    for i in eachindex(Sind)
        sub_idx = Sind[i]
        enz_idx = Eind[i]
        du[sub_idx] = S_subs[i] < 1e-8 && p.q[i] >= 0 ?
            0.0 : p.q[i] * p.E_coli_cellDW * X_tot
        if i == 1
            du[sub_idx] += p.h_release * lysis_term
        end
        mu_avg      = X_tot > 1.0 ?
            (mu_safe * S_cell + mu_eff_l * L_cell) / X_tot : mu_safe
        du[enz_idx] = p.alpha_syn * f[i] * u_cyt[i] -
            (p.beta_deg + mu_avg) * e_enz[i] + 0.001
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

function enforcePositiveDomain!(u)
    for i in eachindex(u)
        u[i] = max(0.0, u[i])
    end
end

function enforceThreshold!(u)
    if u[Sind_S] < CELL_THRESHOLD   u[Sind_S] = 0.0 end
    if u[Sind_I] < CELL_THRESHOLD   u[Sind_I] = 0.0 end
    if u[Sind_L] < CELL_THRESHOLD   u[Sind_L] = 0.0 end
    if u[Pfind]  < PHAGE_THRESHOLD  u[Pfind]  = 0.0 end
end