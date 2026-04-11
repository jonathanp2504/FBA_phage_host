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

    #if haskey(p.fbaModel.reactions, "R_ATPM")
        #p.fbaModel.reactions["R_ATPM"].lower_bound = 0.0
    #end

    # --- 1. GEMEENSCHAPPELIJKE CONSTRAINTS ---
    f = getMonod(S_subs, p)
    R = getRate(f, p)
    v_cyt = getV_cyt(R)
    # MAAK DE VECTOR AAN om limieten op te slaan
    v_max_kin_vector = zeros(length(Sind))
    for i in eachindex(Sind)
        id = p.ex_ids[i]
        conc = S_subs[i]
        if haskey(p.fbaModel.reactions, id)
            e_relatief = e_enz[i] / p.e_max[i]
            v_max_kin = R[i] * e_relatief * v_cyt[i]
            v_max_fysiek = conc / (p.E_coli_cellDW * X_tot * dt_fba + 1e-10)
            # SLA OP voor later gebruik in Scenario B
            v_max_kin_vector[i] = min(v_max_kin, v_max_fysiek)
            # Dit blijft staan zoals je wilde: strikt kinetisch/fysisch
            p.fbaModel.reactions[id].lower_bound = -min(v_max_kin, v_max_fysiek)
        end
    end

    # --- 2. SCENARIO A: NAÏEF ---
    if haskey(p.fbaModel.reactions, p.benz_id)
        p.fbaModel.reactions[p.benz_id].lower_bound = 0.0
        p.fbaModel.reactions[p.benz_id].upper_bound = 0.0
    end
    sol_N = flux_balance_analysis(p.fbaModel, optimizer = HiGHS.Optimizer)
    
    if !isnothing(sol_N) && !isempty(sol_N.fluxes)
        p.mu_N = sol_N.fluxes[p.biomass_id]
        p.q_N  = getFluxes(sol_N, p.ex_ids)
    else 
        p.mu_N = 0.0; p.q_N .= 0.0
    end

    # --- 3. SCENARIO B: LYSOGEEN ---
    #if haskey(p.fbaModel.reactions, p.biomass_id)
        #p.fbaModel.reactions[p.biomass_id].lower_bound = 0.0
    #end
    # Totale import capaciteit in mmol/gDW/h
    totale_import = sum(v_max_kin_vector)
    
    # De eis: de cel MOET een fractie (f_prod) van zijn import omzetten in Benzonase
    geforceerde_benz_flux = totale_import * p.f_prod
    
    if haskey(p.fbaModel.reactions, p.benz_id)
        p.fbaModel.reactions[p.benz_id].lower_bound = 0.0
        p.fbaModel.reactions[p.benz_id].upper_bound = 1000.0
    end

    sol_l = flux_balance_analysis(p.fbaModel, optimizer = HiGHS.Optimizer)
    
    if !isnothing(sol_l) && !isempty(sol_l.fluxes)
        p.mu_l = sol_l.fluxes[p.biomass_id]
        p.q_l  = getFluxes(sol_l, p.ex_ids)
        p.q_benz_l = sol_l.fluxes[p.benz_id]
        println("HOERA: Flux gevonden: ", p.q_benz_l)
    else 
        # Als hij hier komt, is hij nog steeds Infeasible. 
        # Dit betekent dat de cel op dit moment niet genoeg suiker KAN opnemen
        # om die 0.001 te maken. Check je u0 voor de enzymen!
        p.mu_l = 0.0; p.q_l .= 0.0; p.q_benz_l = 0.0
        println("ALARM: Solver faalt zelfs bij 1e-5! Check je aminozuur-balans.")
    end
end