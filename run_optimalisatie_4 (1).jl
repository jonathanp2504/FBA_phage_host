# ============================================================
#  RUN OPTIMALISATIE MODEL 4
#  Voer eerst run_kalibratie_4.jl uit en pas de BENZ_REF en
#  PN_REF waarden aan in de optimizer bestanden voordat je
#  dit bestand runt.
#
#  Optimizer A: optimaliseert MOI + infectietijdstip
#               beginbiomassa vast op FIXED_BIOMASSA (1e9)
#  Optimizer B: optimaliseert infectietijdstip + beginbiomassa
#               MOI vast op FIXED_MOI_B (2.0)
# ============================================================

include("setup_model_4.jl")
include("optimization_weighted_A_model4.jl")
include("optimization_weighted_B_model4.jl")

# --- Optimizer A ---
println("\n=== Optimizer A: MOI + infectietijdstip (vaste biomassa) ===")
resultaat_A = run_optimization_A(p)
println("\nResultaat A:")
println("  MOI              : ", round(resultaat_A.best_moi,     digits=4))
println("  Infectietijdstip : ", round(resultaat_A.best_t_inf,   digits=2), " h")
println("  Beginbiomassa    : ", round(resultaat_A.best_biomass, sigdigits=3), " cellen/L")
println("  Score            : ", round(resultaat_A.best_score,   digits=6))

# --- Optimizer B ---
println("\n=== Optimizer B: infectietijdstip + biomassa (vaste MOI=2.0) ===")
resultaat_B = run_optimization_B(p)
println("\nResultaat B:")
println("  MOI (vast)       : ", resultaat_B.best_moi)
println("  Infectietijdstip : ", round(resultaat_B.best_t_inf,   digits=2), " h")
println("  Beginbiomassa    : ", round(resultaat_B.best_biomass, sigdigits=3), " cellen/L")
println("  Score            : ", round(resultaat_B.best_score,   digits=6))

println("\n=== Optimalisatie Model 4 voltooid ===")
