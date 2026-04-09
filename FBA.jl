include("./parameters.jl")
using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels 
function fbaUpdate!(u, p::Parameters)
    S_subs = u[Sind]
    e_enz  = u[Eind]
    X_tot  = getTotalBiomass(u, p)
    dt_fba = 10/60 

    # --- 1. GEMEENSCHAPPELIJKE CONSTRAINTS ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    v_cyt = getV_cyt(R)

    for i in eachindex(Sind)
        id = p.ex_ids[i]
        conc = S_subs[i]
        if haskey(p.fbaModel.reactions, id)
            e_relatief = e_enz[i] / p.e_max[i]
            v_max_kin = R[i] * e_relatief * v_cyt[i]
            v_max_fysiek = conc / (p.E_coli_cellDW * X_tot * dt_fba + 1e-10)
            
            # Dit blijft staan zoals je wilde: strikt kinetisch/fysisch
            p.fbaModel.reactions[id].lower_bound = -min(v_max_kin, v_max_fysiek)
        end
    end

    # We berekenen nu maar één optimaal groeiscenario op basis van beschikbare suikers
    # Geen Benzonase-eisen meer in de FBA!
    if haskey(p.fbaModel.reactions, p.benz_id)
        p.fbaModel.reactions[p.benz_id].lower_bound = 0.0
        p.fbaModel.reactions[p.benz_id].upper_bound = 0.0
    end

    sol = flux_balance_analysis(p.fbaModel, optimizer = HiGHS.Optimizer)

    if !isnothing(sol) && !isempty(sol.fluxes)
        # De 'ruwe' groei en opname die de cel KAN halen
        p.mu_N = sol.fluxes[p.biomass_id]
        p.q_N  = getFluxes(sol, p.ex_ids)
        
        # Voor de lysogenen gebruiken we dezelfde basis, de 'straf' passen we toe in de ODE
        p.mu_l = p.mu_N 
        p.q_l  = p.q_N
    else 
        p.mu_N = 0.0; p.q_N .= 0.0
        p.mu_l = 0.0; p.q_l .= 0.0
    end
end
