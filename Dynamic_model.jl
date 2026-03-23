using DifferentialEquations
# Include Metabolism struct, which must be available
include("./metabolism.jl") 

mutable struct Bioreactor
    state::ComponentArray{Float64} # gDW/L and mmol/L
    metabolism::Metabolism
    solution::Any
    # Voeg hier eventueel parameters toe voor faag en LamB
end

# De functie die de afgeleiden berekent (de dFBA-stap)
function dynamic_fba_ode(u, p, t)
    # u: Huidige toestand (biomassa, mediumconcentraties)
    # p: Parameters (niet gebruikt in deze basisversie)
    # t: Tijdstip

    # 1. Update de FBA-meting met de huidige mediumconcentraties
    updateMetabolism!(u.metabolism, u.medium)

    # 2. Haal de berekende groeisnelheid en fluxen op
    mu = getGrowthRate(u.metabolism)
    fluxes = getFluxes(u.metabolism)

    # 3. Initialiseer de afgeleiden (diferentiële vergelijkingen)
    du = zero(u) # Maak een ComponentArray met dezelfde structuur

    # 4. Biomassa ODE: dX/dt = mu * X
    du.biomass = mu * u.biomass

    # 5. Medium ODE's: dC_i/dt = - (v_i * X)
    # Dit gaat uit van een yield van 1, wat simpel is voor een kernmodel
    # Voor de opnamefluxen (v_i):
    for key in keys(u.medium)
        # De fluxen zijn al mmol/gDW/h, dus -flux * biomassa geeft dC/dt in mmol/L/h
        du.medium[key] = -fluxes[key] * u.biomass
    end

    return du
end

function simulateBioreactor(b::Bioreactor, duration::Float64)
    # Definieer het ODE-probleem en los het op
    tspan = (0.0, duration)
    prob = ODEProblem(dynamic_fba_ode, b.state, tspan)
    b.solution = solve(prob, Rodas5(), reltol = 1e-3, abstol = 1e-6)
end