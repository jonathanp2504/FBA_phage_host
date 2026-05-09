import OrdinaryDiffEqCore
if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end
using OrdinaryDiffEq
using DelayDiffEq
using UnPack

include("./parameters.jl")
include("./FBA.jl")

# Fysische detectielimiet: onder dit aantal deeltjes/L -> exact nul
const CELL_THRESHOLD  = 100.0   # onder 100 cellen/L -> biologisch irrelevant
const PHAGE_THRESHOLD = 100.0

# ============================================================
#  run() — state vector heeft 13 elementen
#  MOI zit NIET in de toestandsvector.
#  Bereken post-hoc in de main als:
#    MOI = [min(sol.u[i][Pfind] / max(sol.u[i][Sind_S], 1.0), 10.0)
#           for i in eachindex(sol.u)]
# ============================================================
function run(p::Parameters)
    u0         = zeros(13)
    u0[Sind]   = [4.44, 2.337, 5.42, 0.0]
    u0[Eind]   = [0.95, 0.01, 0.01, 0.01]
    u0[Sind_S] = p.startingBiomass

    tspan = (0.0, p.duration)

    infectionCondition(u, t, integrator) = t == p.infection_time
    infectionAffect!(integrator)         = integrator.u[Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints                  = collect(0:1/60:p.duration)
    fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints

    # FBA update + harde celdrempel elke minuut
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
        tstops   = [p.infection_time; fbaUpdateTimepoints],
        callback = CallbackSet(domainCallBack, infectionCallBack, fbaCallBack))
    return solution
end

# ============================================================
#  ODE rechterkant
#
#  Fix 1: Populatievergelijkingen beschermd met celdrempel —
#          cellen onder CELL_THRESHOLD dragen niet bij.
#
#  Fix 2: prob_lys gebaseerd op vertraagde MOI (tau uur geleden)
#          — beslissing wordt genomen bij infectie, niet bij lysis.
#
#  Fix 3: lysis_term enkel actief als er tau uur geleden ook
#          echt vatbare cellen en fagen waren.
# ============================================================
function simulate_dFBA!(du, u, h, p::Parameters, t)::Nothing
    S_subs = u[Sind]
    e_enz  = u[Eind]
    S_cell = u[Sind_S]
    I_cell = u[Sind_I]
    L_cell = u[Sind_L]

    f     = getMonod(S_subs, p)
    R     = getRate(f, p)
    u_cyt = getU_cyt(R, p)

    # Vertraagde toestand tau uur geleden
    u_tau = h(p, t - p.tau)

    # FIX 3: lysis_term alleen als er tau uur geleden echte infectie was
    # (S_tau > drempel en P_tau > drempel, anders geen lysis)
    S_tau = u_tau[Sind_S]
    P_tau = u_tau[Pfind]
    lysis_term = (t > p.tau && S_tau > CELL_THRESHOLD && P_tau > PHAGE_THRESHOLD) ?
        getNewInfectionFlux(u_tau, p) : 0.0

    # FIX 2: prob_lys op basis van MOI tau uur geleden (infectietijdstip)
    phi_beslissing = (t > p.tau && S_tau > CELL_THRESHOLD) ?
        getPhi(u_tau, p) : 0.0
    prob_lys = getProbLys(phi_beslissing)

    X_tot            = getTotalBiomass(u, p)
    totaal_adsorptie = getTotalAdsorptionFlux(u, p)
    nieuwe_infectie  = getNewInfectionFlux(u, p)

    mu_safe  = max(0.0, p.mu)
    mu_eff_l = mu_safe * (1.0 - p.f_prod)

    # --- Benzonase productie ---
    productie_benz = (mu_safe - mu_eff_l) * L_cell * p.Y_benz

    # --- Substraten & enzymen ---
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

    # --- Populatiedynamica ---
    # FIX 1: cellen onder drempel -> geen groei of infectie
    S_actief = S_cell > CELL_THRESHOLD ? S_cell : 0.0
    L_actief = L_cell > CELL_THRESHOLD ? L_cell : 0.0

    du[Sind_S]  = S_actief > 0.0 ? mu_safe * S_actief - nieuwe_infectie : 0.0
    du[Sind_I]  = (1.0 - prob_lys) * nieuwe_infectie - lysis_term
    du[Sind_L]  = mu_eff_l * L_actief + prob_lys * nieuwe_infectie
    du[Pfind]   = p.b * lysis_term - totaal_adsorptie
    du[Benzind] = productie_benz

    return nothing
end

# ============================================================
#  enforcePositiveDomain! — bij negatieve waarden (domainCallback)
# ============================================================
function enforcePositiveDomain!(u)
    for i in eachindex(u)
        u[i] = max(0.0, u[i])
    end
end

# ============================================================
#  enforceThreshold! — wordt elke minuut aangeroepen via fbaAffect!
#  Zet deeltjes onder fysische detectielimiet op exact nul.
#  Dit voorkomt dat de solver blijft integreren op fantoomwaarden.
# ============================================================
function enforceThreshold!(u)
    if u[Sind_S] < CELL_THRESHOLD   u[Sind_S] = 0.0 end
    if u[Sind_I] < CELL_THRESHOLD   u[Sind_I] = 0.0 end
    if u[Sind_L] < CELL_THRESHOLD   u[Sind_L] = 0.0 end
    if u[Pfind]  < PHAGE_THRESHOLD  u[Pfind]  = 0.0 end
end