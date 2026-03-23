
using Plots, SparseArrays
using AbstractFBCModels
using ComponentArrays
gr()

function get_stoich_matrix(model)
    # prefer explicit metabolite list if available
    rxn_ids = collect(keys(model.reactions))
    met_ids = nothing
    if hasproperty(model, :metabolites)
        try
            mlist = getfield(model, :metabolites)
            if isa(mlist, AbstractVector)
                met_ids = [string(x) for x in mlist]
            end
        catch
            met_ids = nothing
        end
    end

    stoich_maps = Vector{Dict{String,Float64}}(undef, length(rxn_ids))
    metset = Set{String}()
    for (j, rid) in enumerate(rxn_ids)
        rxn = model.reactions[rid]
        found = Dict{String,Float64}()
        # common field holding stoichiometry could be a Dict mapping metabolite->coeff
        for fn in fieldnames(typeof(rxn))
            val = nothing
            try
                val = getfield(rxn, fn)
            catch
                continue
            end
            if isa(val, AbstractDict)
                for (k,v) in val
                    found[string(k)] = Float64(v)
                end
                break
            end
        end
        stoich_maps[j] = found
        for k in keys(found)
            push!(metset, k)
        end
    end

    if met_ids === nothing
        met_ids = collect(metset)
    end
    met_index = Dict(m => i for (i,m) in enumerate(met_ids))

    rows = Int[]; cols = Int[]; vals = Float64[]
    for (j, map) in enumerate(stoich_maps)
        for (mid, coeff) in map
            # if metabolite not in met_index, append it
            if !haskey(met_index, mid)
                met_index[mid] = length(met_ids) + 1
                push!(met_ids, mid)
            end
            i = met_index[mid]
            push!(rows, i); push!(cols, j); push!(vals, coeff)
        end
    end
    S = sparse(rows, cols, vals, length(met_ids), length(rxn_ids))
    return S, met_ids, rxn_ids
end

function plot_stoich_heatmap(model; max_met=200, max_rxn=200, top_by_degree=true, sparsity=false, figsize=(900,700))
    S, mets, rxns = get_stoich_matrix(model)
    M, R = size(S)
    # compute degree (number of nonzeros) and sort to bring dense rows/cols first
    rowdeg = mapslices(x->count(!iszero, x), S, dims=2)[:]
    coldeg = mapslices(x->count(!iszero, x), S, dims=1)[:]
    roworder = sortperm(rowdeg, rev=true)
    colorder = sortperm(coldeg, rev=true)

    # select top metabolites/reactions by degree
    imets = roworder[1:min(length(roworder), max_met)]
    irxns = colorder[1:min(length(colorder), max_rxn)]

    if sparsity
        # show sparsity pattern (spy) for quick overview
        p = plot()
        spy(S[imets, irxns]; markersize=1)
        title!("Stoichiometric sparsity pattern (cropped)")
        return p
    end

    A = Array(S[imets, irxns])
    # signed diverging colormap centered at zero
    vmax = maximum(abs.(A))
    if vmax == 0
        vmax = 1.0
    end
    p = heatmap(
        A;
        xlabel = "reactions (ordered)",
        ylabel = "metabolites (ordered)",
        color = :RdBu,
        clims = (-vmax, vmax),
        size = figsize,
        colorbar_title = "stoich coeff",
    )
    return p, mets[imets], rxns[irxns]
end

function flux_time_matrix(b::Bioreactor, times::AbstractVector; reactions_to_track=nothing)
    # compute fluxes at given times; returns matrix R x T and reaction ids
    # reactions_to_track: vector of reaction ids or nothing
    T = length(times)
    # quick first pass to determine candidate reactions if not provided
    if reactions_to_track === nothing
        # use keys of metabolism.fluxes as candidates
        reactions_all = collect(keys(b.metabolism.fluxes))
    else
        reactions_all = collect(reactions_to_track)
    end
    R = length(reactions_all)
    F = zeros(Float64, R, T)
    for (ti,t) in enumerate(times)
        state = b.solution(t)
        mcopy = deepcopy(b.metabolism)
        try
            updateMetabolism!(mcopy, state.medium)
            fluxes = getFluxes(mcopy)
            for (ri, rxn) in enumerate(reactions_all)
                F[ri,ti] = get(fluxes, String(rxn), 0.0)
            end
        catch e
            @warn "FBA failed at time $t: $e"
        end
    end
    return F, reactions_all
end

function plot_flux_heatmap(b::Bioreactor; times=range(0, stop=10, length=80), reactions_to_track=nothing, top_n=40, figsize=(1000,600))
    # sample fewer times by default if long; caller can override
    times_vec = collect(times)
    F, reactions = flux_time_matrix(b, times_vec; reactions_to_track = reactions_to_track)

    # select top_n reactions by maximum absolute flux
    maxabs = maximum(abs.(F), dims=2)[:]
    order = sortperm(maxabs, rev=true)
    sel = order[1:min(top_n, length(order))]
    Fs = F[sel, :]
    rxns_sorted = reactions[sel]

    # center colormap at zero
    vmax = maximum(abs.(Fs))
    if vmax == 0
        vmax = 1.0
    end
    p = heatmap(
        collect(times_vec),
        1:length(rxns_sorted),
        Fs;
        xlabel = "time",
        ylabel = "reaction (top $(length(rxns_sorted)) by |flux|)",
        yticks = (1:length(rxns_sorted), rxns_sorted),
        color = :RdBu,
        clims = (-vmax, vmax),
        size = figsize,
        colorbar_title = "flux (mmol/gDW/h)",
        right_margin = 10Plots.mm,
    )
    return p, rxns_sorted
end
