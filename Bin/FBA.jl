include("./parameters.jl")
using COBREXA
using HiGHS
using AbstractFBCModels
using JuMP
import SBMLFBCModels

struct FbaSolution
    fluxes::Dict{String, Float64}
end

# ============================================================
#  fbaUpdate! — wordt elke tijdsstap aangeroepen vanuit dFBA
#  Benzonase zit NIET in de FBA, productie via f_prod/Y_benz
# ============================================================
function fbaUpdate!(u::Vector{Float64}, p::Parameters)
    vmax::Vector{Float64} = getVmax(u, p)

    sol = solveFba(p.fbaModel, vmax, p)
    if !isnothing(sol) && !isempty(sol.fluxes)
        p.mu = sol.fluxes[p.biomass_id]
        p.q  = getFluxes(sol, p.ex_ids)
    else
        p.mu = 0.0
        p.q .= 0.0
    end
    p.q_benz = 0.0  # niet via FBA berekend
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
            base_name   = reaction_id,
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
        vmax[i]    = R[i] * e_relatief * v_cyt[i]
    end
    return vmax
end

function loadFBAmodel(path)
    model = convert(AbstractFBCModels.CanonicalModel.Model, load_model(path))
    model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0
    return model
end