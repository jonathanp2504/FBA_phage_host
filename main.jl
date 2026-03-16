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
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model("e_coli_core.xml"))

# 2. Cybernetische Parameters 
alpha_syn = 0.2   # Snelheid van enzymsynthese
beta_deg = 0.05   # Snelheid van enzymdegradatie
K_s = 0.005      # Affiniteit uit Tabel 7 (5 mg/L)

include("./parameters.jl")
parameters::Parameters = Parameters(alpha_syn, beta_deg, K_s)

function dFBA_phage_system(du, u, h, p::Parameters, t)
    # INDEX: 
    # 1: X (Biomassa) | 2:4 S (Glc, Mal, Ac) | 5:7 e (Enzymen)
    # 8: S_cell (Vatbaar) | 9: I (Infect) | 10: L (Lysogeen) | 11: P (Fagen)
    
    X = u[1]
    S_subs = u[2:4]
    e_enz  = u[5:7]
    S_cell, I_cell, L_cell, P_phage = u[8:11]
    X_tot = S_cell + I_cell + L_cell

    # --- PARAMETERS UIT TABEL 7 ---
    tau = 0.6          # Latente periode (uren)
    b = 170.0          # Burst size
    alfa_ads = 0.0594  # Adsorptieconstante (L/gDW/h)
    p_pref = [0.8925, 0.08925, 0.01825] # Voorkeurshiërarchie uit thesis
    V_max = [15.0, 13.0, 4.0]           # Opnamesnelheden

    # 1. CYBERNETICA (f, u en v)
    f = [S_subs[i] > 1e-5 ? S_subs[i] / (S_subs[i] + K_s) : 0.0 for i in 1:3]
    R = V_max .* f
    
    # u_cyt bepaalt de synthese van enzymen (e)
    denom = sum(p_pref .* R) + 1e-10
    u_cyt = [(p_pref[i] * R[i]) / denom for i in 1:3]
    
    # v_cyt moduleert de werkelijke opname (v)
    v_cyt = [R[i] / (maximum(R) + 1e-10) for i in 1:3]

    # 2. LYSIS TERMIJN 
    if t < tau
        lysis_term = 0.0
    else
        u_p = h(p, t - tau)
        # Cellen die nu barsten zijn tau geleden geïnfecteerd
        lysis_term = alfa_ads * u_p[8] * u_p[11] 
    end

    # 3. FBA KOPPELING
    for id in keys(model.reactions); if startswith(id, "R_EX_"); model.reactions[id].lower_bound = 0.0; end; end
    
    # De bounds worden beperkt door het huidige enzym-niveau e_enz
    model.reactions["R_EX_glc__D_e"].lower_bound = -R[1] * e_enz[1] * v_cyt[1]
    model.reactions["R_EX_mal__L_e"].lower_bound = -R[2] * e_enz[2] * v_cyt[2]
    model.reactions["R_EX_ac_e"].lower_bound      = -R[3] * e_enz[3] * v_cyt[3]

    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    
    if isnothing(sol) || isnan(sol.fluxes["R_BIOMASS_Ecoli_core_w_GAM"])
        mu, q = 0.0, zeros(3)
    else
        mu = sol.fluxes["R_BIOMASS_Ecoli_core_w_GAM"]
        q = [sol.fluxes["R_EX_glc__D_e"], sol.fluxes["R_EX_mal__L_e"], sol.fluxes["R_EX_ac_e"]]
    end

    # 4. DIFFERENTIAALVERGELIJKINGEN
    du[1] = mu * X_tot 
    
    for i in 1:3
        du[i+1] = q[i] * X_tot 
        # Enzymdynamiek: Synthese (alfa * f * u) - Degradatie (beta * e)
        du[i+4] = p.alpha_syn * f[i] * u_cyt[i] - p.beta_deg * e_enz[i]
    end

    # Faag-Host interactie
    du[8] = mu * S_cell - (alfa_ads * S_cell * P_phage)
    du[9] = (alfa_ads * S_cell * P_phage) - lysis_term
    du[10] = mu * L_cell
    du[11] = (b * lysis_term) - (alfa_ads * S_cell * P_phage)
end

# 5. INITIALISATIE
# [X, Glc, Mal, Ac, e_glc, e_mal, e_ac, S, I, L, P]
u0 = [0.001, 8.0, 8.0, 0.0, 0.97, 0.01, 0.01, 0.001, 0.0, 0.0, 1000.0]
tspan = (0.0, 48.0) # Tijd verlengd naar 48u omdat de start-biomassa lager is

prob = DDEProblem(dFBA_phage_system, u0, (p,t)->u0, tspan, parameters)
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-6)

# 6. UITGEBREID PLOTTEN
# Plot A: Substraatverloop
p1 = plot(sol, idxs=[2,3,4], title="Substraten (g/L)", 
          label=["Glucose" "Maltose" "Acetaat"], lw=2, ylabel="Conc.")

# Plot B: Enzymatische Adaptatie 
p2 = plot(sol, idxs=[5,6,7], title="Enzym-niveaus (Cybernetische e)", 
          label=["e_Glc" "e_Mal" "e_Ac"], lw=2, ls=:dash, ylabel="Relatief niveau")


# Plot C: Host populatie
p3 = plot(sol.t, [u[8]+u[9]+u[10] for u in sol.u], title="Bacteriële Populatie", 
          label="Totaal X", color=:black, lw=3, ylabel="gDW/L")
plot!(p3, sol, idxs=[8,9], label=["Vatbaar (S)" "Infect (I)"], alpha=0.7)


# Plot D: Faag Lambda Dynamiek
p4 = plot(sol, idxs=[11], title="Fagen (P)", yscale=:log10, 
          color=:red, lw=2, label="Phage Lambda", ylabel="Log10(P)")


plot(p1, p2, p3, p4, layout=(4,1), size=(800,1400), margin=5Plots.mm)