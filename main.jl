using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels
using Pkg
#Pkg.add("DifferentialEquations")
#Pkg.add("DelayDiffEq")
using Plots
using DifferentialEquations


# 1. Model laden
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model("./e_coli_core.xml"))

# 2. Cybernetische Parameters 
alpha_syn = 0.2   # Snelheid van enzymsynthese
beta_deg = 0.05   # Snelheid van enzymdegradatie
K_s = [0.0278, 0.0146, 0.0833]      # Affiniteit uit Tabel 7 (mmol/L)
tau = 0.6          # Latente periode (uren)
b = 170.0          # Burst size
alfa_ads = 1e-7  # Adsorptieconstante (L/gDW/h)
p_pref = [0.8925, 0.08925, 0.01825] # Voorkeurshiërarchie uit thesis
V_max = [15.0, 13.0, 4.0]           # Opnamesnelheden (mmol/hr^-1)
E_coli_cellDW = 2.8e-13 # gDW per cel

include("./parameters.jl")
parameters::Parameters = Parameters(alpha_syn, beta_deg, K_s, tau, b, alfa_ads, p_pref, V_max, E_coli_cellDW,7,8,9,10,4:6,1:3)
include("./dFBA_function.jl")
# 5. INITIALISATIE
# [Glc, Mal, Ac, e_glc, e_mal, e_ac, S, I, L, P] (subs in mmol/l)
u0 = [4.44, 2.337, 0.0, 0.97, 0.01, 0.01, 1e6, 0.0, 0.0, 0.0]
tspan = (0.0, 15.0) # Tijd verlengd naar 48u omdat de start-biomassa lager is

infectionCondition(u, t, integrator) = t == 5.0 # h !! VOEG TOE AAN Parameters
infectionAffect!(integrator) = integrator.u[parameters.ind_P] = 1000.0
infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

prob = DDEProblem(dFBA_phage_system, u0, (p,t)->u0, tspan, parameters)
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-6, tstops=[5.0], callback=infectionCallBack )

# 6. UITGEBREID PLOTTEN
# Plot A: Substraatverloop
p1 = plot(sol, idxs=parameters.ind_subs, title="Substraten (mmol/L)", 
          label=["Glucose" "Maltose" "Acetaat"], lw=2, ylabel="Conc.")

# Plot B: Enzymatische Adaptatie 
p2 = plot(sol, idxs=parameters.ind_e, title="Enzym-niveaus (Cybernetische e)", 
          label=["e_Glc" "e_Mal" "e_Ac"], lw=2, ls=:dash, ylabel="Relatief niveau")


# Plot C: Host populatie
p3 = plot(sol.t, [getTotalBiomass(u, parameters) for u in sol.u], title="Bacteriële Populatie", 
          label="Totaal X", color=:black, lw=3, ylabel="gDW/L", ylims=:auto)
plot!(p3, sol, idxs=[parameters.ind_S, parameters.ind_I], label=["Vatbaar (S)" "Infect (I)"], alpha=0.7)


# Plot D: Faag Lambda Dynamiek
p4 = plot(sol, idxs=[parameters.ind_P], title="Fagen (P)", yscale=:log10, 
          color=:red, lw=2, label="Phage Lambda", ylabel="Log10(P)")


plot(p1, p2, p3, p4, layout=(2,2), size=(800,800), margin=5Plots.mm)