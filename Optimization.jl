Pkg.add(["Optimization", "OptimizationOptimJL", "ForwardDiff"])
include("./setup.jl")
using Optimization, OptimizationOptimJL
using DifferentialEquations
include("./Bin/parameters.jl")
include("./Bin/FBA.jl")
include("./Bin/dFBA.jl")

# Zorg dat je hoofdfuncties beschikbaar zijn (pas de filename aan indien nodig)
# include("main_model.jl") 

function run_optimization(p_initial::Parameters, tspan)
    
    # 1. De Objective Functie
    # x bevat de parameters die we willen tunen, bijv:
    # x[1] = MOI (Multiplicity of Infection)
    # x[2] = t_induction (Wanneer voeg je de fagen toe)
    function objective(x, p_opt)
        # Maak een kopie van de parameters en update de waarden
        p_current = remake(p_opt, 
            MOI = x[1], 
            t_induction = x[2]
        )
        
        # Definieer en solve het DDE probleem
        # We gebruiken een try-catch omdat sommige parameter-combinaties 
        # het model kunnen laten crashen (bijv. mu -> 0 te vroeg)
        try
            prob = DDEProblem(dFBA_phage_system, u0, h, tspan, p_current)
            sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-6, abstol=1e-6)
            
            # We willen Benzonase aan het einde MAXIMALISEREN
            # Optimization.jl minimaliseert, dus we nemen de negatieve waarde
            return -sol[Benzind, end] 
        catch
            return 0.0 # Straf voor gecrashte simulaties
        end
    end

    # 2. Setup van de Optimalisatie
    # x0: initiële gok [MOI, t_inductie]
    x0 = [1.0, 3.0] 
    lb = [1e-5, 1.0]  # Ondergrenzen
    ub = [2.0, 15.0] # Bovengrenzen (bijv. max 30 uur)

    opt_func = OptimizationFunction(objective, Optimization.AutoForwardDiff())
    prob_opt = OptimizationProblem(opt_func, x0, p_initial, lb = lb, ub = ub)

    # 3. Solver aanroepen
    # NelderMead is robuust voor dFBA modellen
    sol = solve(prob_opt, NelderMead(), maxiters=100)
    
    return sol
end

# Voorbeeld aanroep:
# result = run_optimization(p, (0.0, 48.0))
# println("Optimale MOI: ", result.u[1], " Optimale tijd: ", result.u[2])