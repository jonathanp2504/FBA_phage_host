using Plots
include("./bioreactor (1).jl")


# Build the combined plot and save a PNG next to this file so you can open it
function plotBioreactor(time::Vector{Float64}, b::Bioreactor)
    pBiomass = plot(time, [b.solution(t).biomass for t in time], xlabel="h", ylabel="gDW/L", labels="biomass")
    pSubstrate = plot(time, [b.solution(t).medium.R_EX_glc__D_e for t in time], xlabel="h", ylabel="mmol/L", labels="glc")
    pSubstrate = plot!(pSubstrate, time, [b.solution(t).medium.R_EX_mal__L_e for t in time], labels="mal")
    p = plot(pBiomass, pSubstrate, size=(1000,600))

    # save output next to this file; non-fatal if it fails
    try
        outpath = joinpath(@__DIR__, "bioreactor_plot.png")
        savefig(p, outpath)
    catch err
        @warn "plotting: failed to save figure" error=err
    end

    return p
end