include("./parameters.jl")
# Plot A: Substraatverloop
function plotSubstrates(sol)
    return plot(sol, idxs=Sind, title="Substraten (mmol/L)", 
          label=["Glucose" "Maltose" "Glycerol" "Acetaat"], lw=2, xlabel= "t [h]", ylabel="Conc.")
end
function plotEnzymes(sol)
    return plot(sol, idxs=Eind, title="Enzym-niveaus (Cybernetische e)", 
          label=["e_Glc" "e_Mal" "e_Glyc" "e_Ac"], lw=2, ls=:dash, xlabel= "t [h]", ylabel="Relatief niveau")
end
function plotPopulation(sol, p)
    X_totaal_data = [getTotalBiomass(u, p) for u in sol.u]

    p3 = plot(sol.t, X_totaal_data, title="Bacteriële Populatie", yscale=:log10, ylims=(1,:auto),
              label="Totaal X", color=:black, lw=3, ylabel="cells [1/L]")
    plot!(p3, sol, idxs=[Nind, Dind, Lind, lind], label=["Naive (N)" "Deciding (D)" "Lytic (I)" "Lysogeen (L)"], alpha=0.7)
    return p3
end
function plotPhages(sol)
    return plot(sol, idxs=[Pfind, Paind], label=["Free phages" "Attached Phages"], title="Fagen (P)", yscale=:log10, ylims=(1,:auto),
        lw=2, xlabel= "t [h]", ylabel="phages [1/L]")
end

function plotAll(sol, p)
    p1 = plotSubstrates(sol)
    p2 = plotEnzymes(sol)
    p3 = plotPopulation(sol, p)
    p4 = plotPhages(sol)
    return plot(p1, p2, p3, p4, 
                  layout=(2, 2), 
                  size=(800, 800), # Iets groter gemaakt voor de leesbaarheid
                  margin=5Plots.mm)
end

#= Plot E: Benzonase Accumulatie
p5 = plot(sol, idxs=Benzind, title="Benzonase Productie", 
          color=:green, lw=2, 
          xlabel="t [h]", ylabel="Conc. [mmol/L]",
          label="Benzonase")

# Plot voor de individuele ladder-treden
p6 = plot(title="Adsorptie Ladder Verdeling", yscale=:log10, 
          xlabel="t [h]", ylabel="cells [1/L]", ylims=(1, :auto))
=#
# Combineer alle plots in een nieuwe layout

#volgende stap is nu productie erin te krijgen en in uw FBA inbouwen door bijvoorbeeld zelf reacties toe te voegen zoals in stoich matrix (bijvoorbeeld met insuline/benzoase (enzym misschien goede eerste stap, moet niet per se aan FBA omdat enzym) --> welke enzymes etc zijn nodig in E. Coli hiervoor en welke substraten moeten voor deze productie worden opgenomen, zie paper in edge)
#branch maken in github ook dat je terugkan als het niet lukt
# als ik dat heb zoek op optimization.jl (zoek package en neem door)