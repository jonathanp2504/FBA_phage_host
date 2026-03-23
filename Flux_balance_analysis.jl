using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels

# 1. Laad het model
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model("e_coli_core.xml"))

# 2. Definieer het medium (Constraints)
# Stel we hebben een mix van Glucose en maltose en lactose beschikbaar, en zuurstof in overvloed
model.reactions["R_EX_glc__D_e"].lower_bound = -20.0   # Glucose opname
model.reactions["R_EX_glc__D_e"].upper_bound = 1000.0
model.reactions["R_EX_mal__L_e"].lower_bound = -20.0   # "Maltose" opname
model.reactions["R_EX_mal__L_e"].upper_bound = 1000.0
model.reactions["R_EX_lac__D_e"].lower_bound = -20.0
model.reactions["R_EX_lac__D_e"].upper_bound = 1000.0
model.reactions["R_EX_o2_e"].lower_bound = -20.0    # Zuurstof
model.reactions["R_EX_o2_e"].upper_bound = 1000.0
model.reactions["R_BIOMASS_Ecoli_core_w_GAM"].upper_bound = 2.0 
#later redenen vinden voor bounds
# 3. Voer de FBA uit
solution = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)

# 4. Resultaten ophalen (veerkrachtig tegen ontbrekende reactienamen)
function get_flux(sol, rxn)
	if haskey(sol.fluxes, rxn)
		return sol.fluxes[rxn]
	else
		@warn "Flux key niet gevonden" reaction=rxn
		return NaN
	end
end

growth = get_flux(solution, "R_BIOMASS_Ecoli_core_w_GAM")
acetate = get_flux(solution, "R_EX_ac_e")
glucose_flux = get_flux(solution, "R_EX_glc__D_e")
maltose_flux = get_flux(solution, "R_EX_mal__L_e")
lac_flux = get_flux(solution, "R_EX_lac__D_e")

println("--- FBA Resultaten ---")
println("Groei snelheid: ", round(growth, digits=4), " 1/h")
println("Glucose flux:   ", round(glucose_flux, digits=2), " mmol/gDW/h")
println("Maltose flux:   ", round(maltose_flux, digits=2), " mmol/gDW/h")
println("Lactose flux:   ", round(lac_flux, digits=2), " mmol/gDW/h")
println("Acetaat flux:   ", round(acetate, digits=2), " mmol/gDW/h")

println("\nBeschikbare reactienamen (eerste 40):")
println(join(sort(collect(keys(solution.fluxes)))[1:40], ", "))