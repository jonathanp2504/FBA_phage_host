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

    # Naïeve cellen: gewone FBA zonder Benzonase
    sol = solveFba(p.fbaModelNaive, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_N = sol.fluxes[p.biomass_id]
        p.q_N  = getFluxes(sol, p.ex_ids)
    else
        p.mu_N = 0.0; p.q_N .= 0.0
    end

    # Lysogene cellen: zelfde FBA als naïef
    # Benzonase kost wordt NIET in de FBA gestopt
    # maar verrekend via f_prod en Y_benz in de ODE
    sol = solveFba(p.fbaModelLysogen, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu_l = sol.fluxes[p.biomass_id]
        p.q_l  = getFluxes(sol, p.ex_ids)
    else
        p.mu_l = 0.0; p.q_l .= 0.0
    end
    p.q_benz_l = 0.0
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
    enzymes    = u[Eind]
    f     = getMonod(substrates, p)
    R     = getRate(f, p)
    v_cyt = getV_cyt(R)
    vmax::Vector{Float64} = zeros(length(substrates))
    for i in eachindex(substrates)
        e_relatief = enzymes[i] / p.e_max[i]
        vmax[i] = R[i] * e_relatief * v_cyt[i]
    end
    return vmax
end

function loadFBAmodel(path)
    model = convert(AbstractFBCModels.CanonicalModel.Model, load_model(path))
    model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0
    return model
end

# addBenzonase! staat hier voor compatibiliteit maar wordt in deze branch
# niet gebruikt omdat Benzonase buiten de FBA wordt berekend
function addBenzonase!(model, benz_stoich)
    model.metabolites["M_benzonase_c"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_benzonase_c"].name = "Benzonase"
    model.metabolites["M_benzonase_c"].compartment = "c"
    model.metabolites["M_benzonase_e"] = AbstractFBCModels.CanonicalModel.Metabolite()
    model.metabolites["M_benzonase_e"].name = "Benzonase (extracellular)"
    model.metabolites["M_benzonase_e"].compartment = "e"
    model.reactions["R_BENZ_prod"] = AbstractFBCModels.CanonicalModel.Reaction()
    model.reactions["R_BENZ_prod"].name = "Benzonase production"
    model.reactions["R_BENZ_prod"].lower_bound = 0.0
    model.reactions["R_BENZ_prod"].upper_bound = 1000.0
    model.reactions["R_BENZ_prod"].stoichiometry = benz_stoich
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