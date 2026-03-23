import Pkg
Pkg.add("DifferentialEquations")
using DifferentialEquations
using COBREXA
using HiGHS
using AbstractFBCModels
using Plots
import SBMLFBCModels
# 1. Laad het model
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model("e_coli_core.xml"))

# 2. Definieer het medium (Constraints)
# Stel we hebben een mix van Glucose en maltose en lactose beschikbaar, en zuurstof in overvloed
model.reactions["R_EX_glc__D_e"].lower_bound = -10.32   # Glucose opname
model.reactions["R_EX_glc__D_e"].upper_bound = 1000.0
model.reactions["R_EX_mal__L_e"].lower_bound = 0   # "Maltose" opname
model.reactions["R_EX_mal__L_e"].upper_bound = 1000.0
model.reactions["R_EX_lac__D_e"].lower_bound = 0
model.reactions["R_EX_lac__D_e"].upper_bound = 1000.0
model.reactions["R_EX_o2_e"].lower_bound = -20.0    # Zuurstof
model.reactions["R_EX_o2_e"].upper_bound = 1000.0
model.reactions["R_BIOMASS_Ecoli_core_w_GAM"].upper_bound = 2.0 
#later redenen vinden voor bounds #beter in parameter p van functie uw model meegeven (niet in functie zelf definieren)
function f(u,p,t)
    println(u)
    du = zeros(2)
    model.reactions["R_EX_glc__D_e"].lower_bound = -10.32* (u[2])/((u[2] + 0.02775))   # Glucose opname afhankelijk van huidige concentratie
    println(model.reactions["R_EX_glc__D_e"].lower_bound)
    solution = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    if typeof(solution) == Nothing 
        qs=0.0
        mu = 0.0
    else
        qs = solution.fluxes["R_EX_glc__D_e"] # specifieke glucose opname
        mu = solution.fluxes["R_BIOMASS_Ecoli_core_w_GAM"]
    end
    du[1] = u[1] * mu 
    du[2] = (u[1]*0.3e-12) * qs 
    return du
end
function iod(u, p, t)
    return any(x -> x < 0, u)
end
u0= [1e10, 10]
tspan = (0.0, 10.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob,Tsit5(), isoutofdomain=iod)
using Plots; plot(sol,plotdensity= 1000, idxs = [(0,1)], labels=["biomass", "substrate"], xlabel="time", ylabel="concentration")
using Plots; plot(sol,plotdensity= 1000, idxs = [(0,2)], labels=["biomass", "substrate"], xlabel="time", ylabel="concentration")

# 1. Laad het model
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model("e_coli_core.xml"))

# 2. Definieer het medium (Constraints)
# Stel we hebben een mix van Glucose en maltose en lactose beschikbaar, en zuurstof in overvloed
model.reactions["R_EX_glc__D_e"].lower_bound = 0   # Glucose opname
model.reactions["R_EX_glc__D_e"].upper_bound = 1000.0
model.reactions["R_EX_mal__L_e"].lower_bound = 0   # "Maltose" opname
model.reactions["R_EX_mal__L_e"].upper_bound = 1000.0
model.reactions["R_EX_lac__D_e"].lower_bound = 0
model.reactions["R_EX_lac__D_e"].upper_bound = 1000.0
model.reactions["R_EX_o2_e"].lower_bound = -20.0    # Zuurstof
model.reactions["R_EX_o2_e"].upper_bound = 1000.0
model.reactions["R_BIOMASS_Ecoli_core_w_GAM"].upper_bound = 2.0 
#later redenen vinden voor bounds
# 3. Voer de FBA uit
solution = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)

# 4. Resultaten ophalen
growth = solution.fluxes["R_BIOMASS_Ecoli_core_w_GAM"]
acetate = solution.fluxes["R_EX_ac_e"]
glucose_flux = solution.fluxes["R_EX_glc__D_e"]
maltose_flux = solution.fluxes["R_EX_mal__L_e"]
lac_flux = solution.fluxes["R_EX_lac__D_e"]

println("--- FBA Resultaten ---")
println("Groei snelheid: ", round(growth, digits=4), " h⁻¹")
println("Glucose flux:   ", round(glucose_flux, digits=2), " mmol/gDW/h")
println("Maltose flux:   ", round(maltose_flux, digits=2), " mmol/gDW/h")
println("Lactose flux:   ", round(lac_flux, digits=2), " mmol/gDW/h")
println("Acetaat flux:   ", round(acetate, digits=2), " mmol/gDW/h")