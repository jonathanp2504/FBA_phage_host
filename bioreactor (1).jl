using DifferentialEquations
using ComponentArrays
include("./metabolism.jl")
mutable struct Bioreactor
    u0::ComponentArray{Float64}
    metabolism::Metabolism

    solution

    function Bioreactor(u0::ComponentArray{Float64}, metabolism::Metabolism)
        bioreactor = new(
            u0, metabolism
        )
        return bioreactor
    end
end

function simulateBioreactor(bioreactor::Bioreactor, duration::Number)
    iod(u, p, t) = any(u.<0) # iod = is out of domain --> stops simulation if condition is met
    # setup DDE problem
    bioreactor.solution = solve(ODEProblem(bioreactorODEFunction, 
        bioreactor.u0,
        (0,duration), bioreactor; alg=Tsit5()), progress = true, isoutofdomain = iod, dtmax=1e-1)
end


# https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dde_example/
# This is the ODE Function!!!
# each step 
# 1) update the FBA model (metabolism) using the current medium 
# 2) calculate the rate of change for the biomass and medium (the solver uses du to calculate the state of the next step) 
function bioreactorODEFunction(du, u, b::Bioreactor, t)
    updateMetabolism!(b.metabolism, u.medium) 
    du.biomass = getGrowthRate(b.metabolism) * u.biomass
    du.medium = getFluxes(b.metabolism) .* u.biomass
    nothing
end