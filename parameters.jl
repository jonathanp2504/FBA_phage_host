struct Parameters 
    alpha_syn::Float64 
    beta_deg::Float64
    K_s::Vector{Float64}
    tau::Float64
    b::Float64
    alfa_ads::Float64
    p_pref::Vector{Float64}
    V_max::Vector{Float64}
    E_coli_cellDW::Float64

    ind_S::Int64
    ind_I::Int64
    ind_L::Int64
    ind_P::Int64
    ind_e::Vector{Int64}
    ind_subs::Vector{Int64}
end

function getTotalBiomass(u, p::Parameters)
    return u[p.ind_S] + u[p.ind_I] + u[p.ind_L]
end