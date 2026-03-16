using DifferentialEquations
using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels
using Plots

#  Model laden en converteren
raw_model = load_model("iJO1366.xml")
model = convert(AbstractFBCModels.CanonicalModel.Model, raw_model)

# Identificeer alle exchanges
all_exchanges = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]

# Cybernetische Parameters (Tran)
alfa, beta, K_e, K_s = 1.5, 1.2, 0.05, 0.005

# Molmassa's voor de koppeling (g/mol)
MW = [180.16, 342.3, 92.09, 60.05] # Glc, Mal, Glyc, Ac

function dFBA_iml_cybernetic_system(du, u, p, t)
    X, S, e = u[1], u[2:5], u[6:9]
    
    biomass_id = "R_BIOMASS_Ec_iJO1366_core_53p95M"
    ex_ids = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]

    # --- 3. MODEL SETUP ---
    # Reset
    for id in all_exchanges; model.reactions[id].lower_bound = 0.0; end
    
    # UITGEBREIDE MINERALEN & COFACTOREN
    essentials = [
        "R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e", 
        "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e",
        "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e", 
        "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e",
        "R_EX_thi_e", "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"
    ]
    
    for ess in essentials
        if haskey(model.reactions, ess)
            model.reactions[ess].lower_bound = -100.0
        end
    end
    
    # Sta uitstoot toe
    for id in all_exchanges; model.reactions[id].upper_bound = 1000.0; end

    # CYBERNETICA 
    V_max = [15.0, 13.0, 11.0, 4.0] 
    p_pref = [0.8925, 0.08925, 0.009125, 0.009125] 

    f = [S[i] > 1e-4 ? S[i] / (S[i] + K_s) : 0.0 for i in 1:4]
    R = V_max .* f
    
    denom = sum(p_pref .* R) + 1e-10
    u_cyt = [(p_pref[i] * R[i]) / denom for i in 1:4]
    v_cyt = [R[i] / (maximum(R) + 1e-10) for i in 1:4]

    #  FBA 
    for i in 1:4
        if haskey(model.reactions, ex_ids[i])
            model.reactions[ex_ids[i]].lower_bound = -R[i] * e[i] * v_cyt[i]
        end
    end

    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    
    # De structuur die je wilde behouden
    if isnothing(sol) || !haskey(sol.fluxes, biomass_id)
        mu = 0.0
        q = zeros(4)
    else
        mu = sol.fluxes[biomass_id]
        q = [haskey(sol.fluxes, id) ? sol.fluxes[id] : 0.0 for id in ex_ids]
    end

    # ODE DEFINITIES
    du[1] = mu * X 
    for i in 1:4
        du[i+1] = q[i] * X * (MW[i] / 1000.0) 
        du[i+5] = alfa * f[i] * u_cyt[i] - beta * e[i]
    end
end

# INITIALISATIE & SOLVE
u0 = [0.01, 8.0, 8.0, 5.0, 0.0, 0.95, 0.01, 0.01, 0.01] 
tspan = (0.0, 48.0)

prob = ODEProblem(dFBA_iml_cybernetic_system, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-3, abstol=1e-3)

# 8. PLOTTEN
p1 = plot(sol, idxs=[2,3,4,5], label=["Glc" "Mal" "Glyc" "Ac"], title="Substraten (g/L)", lw=2)
p2 = plot(sol, idxs=[6,7,8,9], label=["e_Glc" "e_Mal" "e_Glyc" "e_Ac"], title="Enzymen (e)", ls=:dash)
p3 = plot(sol, idxs=[1], label="Biomassa", title="Biomassa (gDW/L)", color=:black, lw=2, ylims=:auto)

l = @layout [a; b; c]
plot(p1, p2, p3, layout=l, size=(800,1000), margin=5Plots.mm)

