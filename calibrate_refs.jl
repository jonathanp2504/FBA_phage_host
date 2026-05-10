# ============================================================
#  KALIBRATIESCRIPT voor BENZ_REF en PN_REF
#
#  Draai dit script EENMALIG vóór je optimizer om de
#  referentiewaarden in te stellen.
#
#  Plak de uitvoer in optimization_weighted_A.jl en
#  optimization_weighted_B.jl als BENZ_REF en PN_REF.
# ============================================================

println("=== Kalibratie van referentiewaarden ===\n")

sol = run(p)

if SciMLBase.successful_retcode(sol)

    # --- Benzonase ---
    max_benz = maximum(sol[Benzind, :])
    println("Maximale Benzonase tijdens simulatie : $(max_benz) mmol/L")
    println("Stel BENZ_REF in op               : $(round(max_benz, sigdigits=1))")

    # --- P/N verhouding over de tijd ---
    P_tijdsreeks  = [sol.u[i][Pfind] for i in eachindex(sol.u)]
    N_tijdsreeks  = [sol.u[i][Nind]  for i in eachindex(sol.u)]
    PN_tijdsreeks = [N_tijdsreeks[i] > 1.0 ? P_tijdsreeks[i]/N_tijdsreeks[i] : 0.0
                     for i in eachindex(sol.t)]

    max_PN = maximum(PN_tijdsreeks)
    end_PN = PN_tijdsreeks[end]

    println("\nMaximale P/N verhouding tijdens simulatie : $(round(max_PN, sigdigits=2))")
    println("Eindwaarde P/N verhouding               : $(round(end_PN, sigdigits=2))")

    # Gebruik maximale P/N als referentie om deling door nul te vermijden
    # wanneer alle cellen dood zijn op het einde van de simulatie
    PN_ref = max_PN > 0.0 ? max_PN : 1.0
    println("Stel PN_REF in op (maximum)             : $(round(PN_ref, sigdigits=1))")

    println("""
\n=== Kopieer dit naar je optimizer ===
const BENZ_REF_A  = $(round(max_benz, sigdigits=1))
const PN_REF_A    = $(round(PN_ref,   sigdigits=1))
const BENZ_REF_B  = $(round(max_benz, sigdigits=1))
const PN_REF_B    = $(round(PN_ref,   sigdigits=1))
""")

    # --- Verificatieplot ---
    benz_norm_ts = [sol.u[i][Benzind] / max(max_benz, 1e-10)
                    for i in eachindex(sol.t)]
    PN_norm_ts   = [PN_tijdsreeks[i]  / max(PN_ref, 1e-10)
                    for i in eachindex(sol.t)]
    obj_ts       = 0.8 .* benz_norm_ts .- 0.05 .* PN_norm_ts

    p_check = plot(sol.t, [benz_norm_ts PN_norm_ts],
        label     = ["Benzonase (genorm.)" "P/N (genorm.)"],
        lw        = 2,
        color     = [:green :red],
        xlabel    = "t [h]",
        ylabel    = "Genormaliseerde waarde",
        title     = "Kalibratie: beide termen op dezelfde schaal",
        legend    = :topright)
    plot!(p_check, sol.t, obj_ts,
        label     = "Gewogen objectief (0.8B - 0.05PN)",
        lw        = 2,
        color     = :blue,
        linestyle = :dash)

    display(p_check)
    savefig(p_check, "kalibratie_referentiewaarden.png")
    println("Kalibratiefiguur opgeslagen: kalibratie_referentiewaarden.png")

else
    println("Simulatie gefaald — controleer je parameters")
end