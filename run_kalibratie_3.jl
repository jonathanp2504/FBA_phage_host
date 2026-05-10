# ============================================================
#  RUN KALIBRATIE MODEL 3
#  Laadt model 3 en berekent de referentiewaarden voor de
#  optimizer. Kopieer de uitvoer naar:
#    optimization_weighted_A.jl -> BENZ_REF_A en PN_REF_A
#    optimization_weighted_B.jl -> BENZ_REF_B en PN_REF_B
# ============================================================

include("setup_model_3.jl")
include("calibrate_refs.jl")
