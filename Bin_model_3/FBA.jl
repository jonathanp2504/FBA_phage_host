include("./parameters.jl")
using COBREXA
using HiGHS
using AbstractFBCModels
using JuMP
import SBMLFBCModels 

struct FbaSolution
    fluxes::Dict{String, Float64}
end

function fbaUpdate!(u::Vector{Float64}, p::Parameters)
    vmax::Vector{Float64} = getVmax(u, p)
    ## Solve FBA for naive cells
    sol = solveFba(p.fbaModelNaive, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_N = sol.fluxes[p.biomass_id]
        p.q_N  = getFluxes(sol, p.ex_ids)
    else 
        p.mu_N = 0.0; p.q_N .= 0.0
    end
    ## Solve FBA for lysogenic cells
    # Totale import capaciteit in mmol/gDW/h
    totale_import = sum(vmax)
    
    # De eis: de cel MOET een fractie (f_prod) van zijn import omzetten in Benzonase
    geforceerde_benz_flux = totale_import * p.f_prod
    setLowerBound!(p.fbaModelLysogen, p.benz_id, geforceerde_benz_flux)

    sol = solveFba(p.fbaModelLysogen, vmax, p)
    
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_l = sol.fluxes[p.biomass_id]
        p.q_l  = getFluxes(sol, p.ex_ids)
        p.q_benz_l = sol.fluxes[p.benz_id]
        #println("HOERA: Flux gevonden: ", p.q_benz_l)
    else 
        # Als hij hier komt, is hij nog steeds Infeasible. 
        # Dit betekent dat de cel op dit moment niet genoeg suiker KAN opnemen
        # om die 0.001 te maken. Check je u0 voor de enzymen!
        p.mu_l = 0.0; p.q_l .= 0.0; p.q_benz_l = 0.0
    end
end

function buildFbaCache(model, exchange_ids::Vector{String}, biomass_id::String; benz_id::Union{Nothing, String}=nothing)
    optimizer = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(optimizer)

    reaction_vars = Dict{String, JuMP.VariableRef}()
    for (reaction_id, reaction) in model.reactions
        reaction_vars[reaction_id] = JuMP.@variable(
            optimizer,
            lower_bound = reaction.lower_bound,
            upper_bound = reaction.upper_bound,
            base_name = reaction_id,
        )
    end

    balances = Dict(metabolite_id => JuMP.AffExpr() for metabolite_id in keys(model.metabolites))
    for (reaction_id, reaction) in model.reactions
        reaction_var = reaction_vars[reaction_id]
        for (metabolite_id, coefficient) in reaction.stoichiometry
            JuMP.add_to_expression!(balances[metabolite_id], coefficient, reaction_var)
        end
    end

    for balance in values(balances)
        JuMP.@constraint(optimizer, balance == 0.0)
    end

    JuMP.@objective(optimizer, Max, reaction_vars[biomass_id])

    tracked_ids = copy(exchange_ids)
    push!(tracked_ids, biomass_id)
    if !isnothing(benz_id)
        push!(tracked_ids, benz_id)
    end
    unique!(tracked_ids)

    return FbaCache(optimizer, reaction_vars, copy(exchange_ids), tracked_ids, biomass_id, benz_id)
end

function setLowerBound!(cache::FbaCache, reaction_id::String, value::Float64)
    JuMP.set_lower_bound(cache.reaction_vars[reaction_id], value)
end

function setUpperBound!(cache::FbaCache, reaction_id::String, value::Float64)
    JuMP.set_upper_bound(cache.reaction_vars[reaction_id], value)
end

function solveFba(cache::FbaCache, vmax::Vector{Float64}, p::Parameters)
    # Update alleen de variabele bounds op de reeds gebouwde LP.
    for i in eachindex(vmax)
        setLowerBound!(cache, cache.exchange_ids[i], -vmax[i])
    end

    JuMP.optimize!(cache.optimizer)
    if !JuMP.is_solved_and_feasible(cache.optimizer; dual=false)
        return nothing
    end

    fluxes = Dict{String, Float64}()
    for reaction_id in cache.tracked_ids
        fluxes[reaction_id] = JuMP.value(cache.reaction_vars[reaction_id])
    end

    return FbaSolution(fluxes)
end

function getVmax(u, p::Parameters)
    substrates = u[Sind]
    enzymes = u[Eind]
    f = getMonod(substrates, p)
    R = getRate(f, p)
    v_cyt = getV_cyt(R)
    vmax::Vector{Float64} = zeros(length(substrates))
    for i in eachindex(substrates)
        e_relatief = enzymes[i] / p.e_max[i]
        vmax[i] = R[i] * e_relatief * v_cyt[i]
    end
    return vmax
end

##########################################
# Loading and modifying the FBA model
##########################################
function loadFBAmodel(path)
    model = convert(AbstractFBCModels.CanonicalModel.Model,  load_model(path))
    model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0 # bacterial growth is constrained!!!
    return model
end

function addBenzonase!(model, benz_stoich)
    # Voeg metabolieten toe
    model.metabolites["M_benzonase_c"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_benzonase_c"].name = "Benzonase"
    model.metabolites["M_benzonase_c"].compartment = "c"

    model.metabolites["M_benzonase_e"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_benzonase_e"].name = "Benzonase (extracellular)"
    model.metabolites["M_benzonase_e"].compartment = "e"

    # Voeg productie reactie toe
    model.reactions["R_BENZ_prod"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_BENZ_prod"].name = "Benzonase production"
    model.reactions["R_BENZ_prod"].lower_bound = 0.0
    model.reactions["R_BENZ_prod"].upper_bound = 1000.0
    model.reactions["R_BENZ_prod"].stoichiometry = benz_stoich

    # Voeg export en exchange toe
    model.reactions["R_BENZ_export"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_BENZ_export"].stoichiometry = Dict("M_benzonase_c" => -1.0, "M_benzonase_e" => 1.0)
    model.reactions["R_BENZ_export"].lower_bound = 0.0
    model.reactions["R_BENZ_export"].upper_bound = 1000.0

    model.reactions["R_EX_benz_e"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_EX_benz_e"].stoichiometry = Dict("M_benzonase_e" => -1.0)
    model.reactions["R_EX_benz_e"].lower_bound = 0.0
    model.reactions["R_EX_benz_e"].upper_bound = 1000.0
    return model
end


