using Printf

# Load project code
include("bioreactor (1).jl")

# Build same scenario as main
metabolism = Metabolism("e_coli_core.xml")
medium = ComponentArray{Float64}(R_EX_glc__D_e=500, R_EX_mal__L_e=250)
u = ComponentArray{Float64}(biomass=1e-6, medium = medium)
b = Bioreactor(u, metabolism)

duration = 10
simulateBioreactor(b, duration)

# sample solution
time = collect(0:0.01:duration)
biom = [b.solution(t).biomass for t in time]
glc = [b.solution(t).medium.R_EX_glc__D_e for t in time]
mal = [b.solution(t).medium.R_EX_mal__L_e for t in time]

# basic stats
init_b = biom[1]
final_b = biom[end]
fold = final_b / init_b

# estimate exponential growth rate by linear fit on log(biomass)
using Statistics
valid = biom .> 0
tvec = time[valid]
logb = log.(biom[valid])
# fit slope m in logb = m*t + c
m = sum((tvec .- mean(tvec)).*(logb .- mean(logb))) / sum((tvec .- mean(tvec)).^2)
# doubling time (if m>0) = ln(2)/m
doubling = m > 0 ? log(2)/m : NaN

# substrate changes
init_glc = glc[1]
final_glc = glc[end]
init_mal = mal[1]
final_mal = mal[end]

# time when glucose drops below 50%
half_glc = init_glc * 0.5
idx_half = findfirst(x -> x <= half_glc, glc)
thalf = isnothing(idx_half) ? NaN : time[idx_half]

@printf("Biomass initial: %.5e, final: %.5e, fold change: %.2f\n", init_b, final_b, fold)
@printf("Estimated growth rate m = %.4f 1/h, doubling time = %.2f h\n", m, doubling)
@printf("Glucose initial: %.3f mmol/L, final: %.3f mmol/L\n", init_glc, final_glc)
@printf("Malate initial: %.3f mmol/L, final: %.3f mmol/L\n", init_mal, final_mal)
if !isnan(thalf)
    @printf("Glucose reaches half initial (~%.3f) at t = %.2f h\n", half_glc, thalf)
else
    @printf("Glucose did not drop to half initial within %.1f h\n", duration)
end

# Print a small selected table (every 1 hour)
@printf("\nTime\tBiomass\tGlc\tMal\n")
for (i,t) in enumerate(time)
    if isapprox(t, round(t); atol=1e-6)
        @printf("%.0f\t%.5e\t%.3f\t%.3f\n", t, biom[i], glc[i], mal[i])
    end
end
