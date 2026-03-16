include("./bioreactor (1).jl")

# build metabolism 
metabolism::Metabolism = Metabolism("e_coli_core.xml")

# build medium 
medium = ComponentArray{Float64}(R_EX_glc__D_e=500, R_EX_mal__L_e=250)

# Initial state
u = ComponentArray{Float64}(biomass=1e-6, medium = medium) # gDW/L ; mmol/L
b::Bioreactor = Bioreactor(u, metabolism);

duration = 10
simulateBioreactor(b, duration)

include("./plotting.jl")
# capture the plot object (plotting.jl will also save a PNG next to itself)
p = plotBioreactor(collect(0:0.01:duration), b)
display(p)




