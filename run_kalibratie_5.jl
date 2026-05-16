# ============================================================
#  RUN KALIBRATIE MODEL 5
#  Laadt model 5 en berekent de referentiewaarden voor de
#  optimizer. Kopieer de uitvoer naar:
#    optimization_weighted_A.jl -> BENZ_REF_A en PN_REF_A
#    optimization_weighted_B.jl -> BENZ_REF_B en PN_REF_B
# ============================================================

include("setup_model_5.jl")
include("Model5.jl")
include("calibrate_refs.jl")
lys_max  = maximum([sol.u[i][Model5.lind]   for i in eachindex(sol.u)])
benz_max = maximum([sol.u[i][Model5.Benzind] for i in eachindex(sol.u)])
pn_max   = maximum([sol.u[i][Model5.Pfind] / max(sol.u[i][Model5.Nind], 1.0)
                    for i in eachindex(sol.u)])

println("Max lysogenen : ", lys_max)
println("Max Benzonase : ", benz_max)
println("Max P/N       : ", pn_max)

