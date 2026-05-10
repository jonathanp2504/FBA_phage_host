# ============================================================
#  HOOFDSTUK 4: Burst size — vast, lineair of metabolisch?
#  Model 3 (vast b) vs Model 4 (lineair b(µ)) vs Model 5 (FBA b)
#
#  Analyses:
#   4a. Burst size tijdsverloop: Model 3 vs 4 vs 5
#   4b. Benzonase productie: vergelijking drie modellen
#   4c. Faagpopulatiedynamica: effect van variabele burst size
#   4d. Optimizer score vergelijking
# ============================================================

using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels, OrdinaryDiffEqCore

if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end

# ============================================================
#  Modellen laden
#  LET OP: elk model heeft zijn eigen Parameters struct
#  Gebruik aparte Julia sessies of modules voor productie
#  Hier worden ze sequentieel geladen en gesimuleerd
# ============================================================

model_path    = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn     = 2.0;  beta_deg = 0.5
K_s           = [0.0278, 0.0146, 0.0543, 0.0833]
tau           = 1.0
p_pref        = [0.8925, 0.08925, 0.008925, 0.008925]
V_max         = [0.0, 3.75, 0.0, 4.0]
E_coli_cellDW = 1e-12
exchange_ids  = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values     = [180.16, 342.3, 92.09, 60.05]
h_release     = 6.0e-12;  duration = 20.0
mu_max_vector = [1.33, 1.26, 1.10, 0.29]
e_max_vector  = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
N0=1e9; t_inf=2.0; moi=2.0; f_prod=0.0015

# Gemeenschappelijke FBA modellen
naiveModel   = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
all_ex_ids   = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
naiveFba     = buildFbaCache(naiveModel,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba   = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                              benz_id="R_BENZ_prod")

# ============================================================
#  Simuleer Model 3 (vast b=170)
# ============================================================
println("=== Simulatie Model 3 (vast b=170) ===")
include("./Bin_model_3/parameters.jl")
include("./Bin_model_3/FBA.jl")
include("./Bin_model_3/dFBA.jl")

p3 = Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, 170.0, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba, lysogenFba, 0.0, zeros(4), 0.0, zeros(4), 0.0,
    1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
    "R_BENZ_prod", 0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)
sol3 = run(p3)

t3   = sol3.t
B3   = [sol3.u[i][Benzind] for i in eachindex(sol3.u)]
Pf3  = [sol3.u[i][Pfind]  for i in eachindex(sol3.u)]
Lys3 = [sol3.u[i][lind]   for i in eachindex(sol3.u)]
N3   = [sol3.u[i][Nind]   for i in eachindex(sol3.u)]
b_tijds_m3 = fill(170.0, length(t3))   # vaste burst size
println("  Max Benzonase: $(round(maximum(B3), sigdigits=3)) mmol/L")

# ============================================================
#  Simuleer Model 4 (lineaire b(µ) = b0 + kb*µN)
# ============================================================
println("=== Simulatie Model 4 (lineaire burst size) ===")
include("./Bin_model_4/parameters.jl")
include("./Bin_model_4/FBA.jl")
include("./Bin_model_4/dFBA.jl")

b0 = 50.0; kb = 95.2
p4 = Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, b0, kb, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba, lysogenFba, 0.0, zeros(4), 0.0, zeros(4), 0.0,
    1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
    "R_BENZ_prod", 0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)
sol4 = run(p4)

t4   = sol4.t
B4   = [sol4.u[i][Benzind] for i in eachindex(sol4.u)]
Pf4  = [sol4.u[i][Pfind]  for i in eachindex(sol4.u)]
Lys4 = [sol4.u[i][lind]   for i in eachindex(sol4.u)]
N4   = [sol4.u[i][Nind]   for i in eachindex(sol4.u)]
# Burst size schatting: b0 + kb * mu_N (benadering via interpolatie)
b_tijds_m4 = b0 .+ kb .* 0.8 .* ones(length(t4))  # aanname mu_N ≈ 0.8 gemiddeld
println("  Max Benzonase: $(round(maximum(B4), sigdigits=3)) mmol/L")

# ============================================================
#  Simuleer Model 5 (FBA burst size)
# ============================================================
println("=== Simulatie Model 5 (FBA burst size) ===")
include("./Bin_model_5/parameters.jl")
include("./Bin_model_5/FBA.jl")
include("./Bin_model_5/dFBA.jl")

lyticModel = addPhage!(loadFBAmodel(model_path), phage_stoich)
lyticFba   = buildLyticFbaCache(lyticModel, exchange_ids,
                "R_BIOMASS_Ec_iJO1366_core_53p95M", "R_PHAGE_prod")

p5 = Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba, lysogenFba, lyticFba,
    0.0, zeros(4), 0.0, zeros(4), 0.0, 0.0,
    1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
    "R_BENZ_prod", "R_PHAGE_prod",
    0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)
sol5 = run(p5)

t5   = sol5.t
B5   = [sol5.u[i][Benzind] for i in eachindex(sol5.u)]
Pf5  = [sol5.u[i][Pfind]  for i in eachindex(sol5.u)]
Lys5 = [sol5.u[i][lind]   for i in eachindex(sol5.u)]
N5   = [sol5.u[i][Nind]   for i in eachindex(sol5.u)]
println("  Max Benzonase: $(round(maximum(B5), sigdigits=3)) mmol/L")

# ============================================================
#  4a. Burst size tijdsverloop
# ============================================================
println("\n=== 4a: Burst size tijdsverloop ===")

fig4a = plot(t3, b_tijds_m3,
    label="Model 3: vast b=170", color=:steelblue, lw=2, linestyle=:dash,
    xlabel="t [h]", ylabel="Burst size b [fagen/cel]",
    title="4a: Burst size over de tijd",
    ylims=(0, 250), legend=:topright)
plot!(fig4a, t4, b_tijds_m4,
    label="Model 4: lineair b(µ)", color=:darkorange, lw=2)
annotate!(fig4a, duration*0.7, 180,
    text("Model 5: b uit FBA\n(zie tekst)", 8, :darkgreen))
savefig(fig4a, "h4a_burst_size_tijdsverloop.png")
println("  Figuur opgeslagen: h4a_burst_size_tijdsverloop.png")

# ============================================================
#  4b. Benzonase productie: drie modellen
# ============================================================
println("\n=== 4b: Benzonase vergelijking ===")

fig4b = plot(t3, B3, label="Model 3 (vast b)", color=:steelblue, lw=2)
plot!(fig4b, t4, B4, label="Model 4 (lineair b)", color=:darkorange, lw=2,
    linestyle=:dash)
plot!(fig4b, t5, B5, label="Model 5 (FBA b)", color=:darkgreen, lw=2,
    linestyle=:dot)
xlabel!(fig4b, "t [h]"); ylabel!(fig4b, "Benzonase [mmol/L]")
title!(fig4b, "4b: Benzonase productie — Model 3 vs 4 vs 5")
savefig(fig4b, "h4b_benzonase_drie_modellen.png")
println("  Figuur opgeslagen: h4b_benzonase_drie_modellen.png")

# ============================================================
#  4c. Faagpopulatiedynamica
# ============================================================
println("\n=== 4c: Faagpopulatiedynamica ===")

fig4c = plot(t3, Pf3, label="Model 3 (vast b)", color=:steelblue, lw=2,
    yscale=:log10, ylims=(1,:auto),
    xlabel="t [h]", ylabel="Vrije fagen/L",
    title="4c: Faagpopulatie — drie modellen")
plot!(fig4c, t4, Pf4, label="Model 4 (lineair b)",
    color=:darkorange, lw=2, linestyle=:dash)
plot!(fig4c, t5, Pf5, label="Model 5 (FBA b)",
    color=:darkgreen, lw=2, linestyle=:dot)
savefig(fig4c, "h4c_fagen_drie_modellen.png")
println("  Figuur opgeslagen: h4c_fagen_drie_modellen.png")

# ============================================================
#  4d. Samenvattingstabel: kernresultaten per model
# ============================================================
println("\n=== 4d: Samenvattingstabel ===")

function pn_eind(sol, Pfind_, Nind_)
    P_f = sol.u[end][Pfind_]; N_f = sol.u[end][Nind_]
    return N_f > 1.0 ? P_f/N_f : Inf
end

benz_ref = max(maximum(B3), maximum(B4), maximum(B5))
pn_ref   = max(pn_eind(sol3, Pfind, Nind),
               pn_eind(sol4, Pfind, Nind),
               pn_eind(sol5, Pfind, Nind))

results = [
    ("Model 3", maximum(B3), pn_eind(sol3, Pfind, Nind)),
    ("Model 4", maximum(B4), pn_eind(sol4, Pfind, Nind)),
    ("Model 5", maximum(B5), pn_eind(sol5, Pfind, Nind)),
]

println("\n  Model   | Max Benz [mmol/L] | P/N eind  | Score (0.8B-0.2PN)")
println("  --------|-------------------|-----------|--------------------")
for (naam, benz, pn) in results
    score = 0.8*(benz/benz_ref) - 0.2*(pn/pn_ref)
    @printf("  %-7s | %-17.5f | %-9.2e | %.4f\n", naam, benz, pn, score)
end

# Barplot samenvatting
fig4d = groupedbar(
    ["Model 3" "Model 4" "Model 5"],
    hcat([r[2] for r in results]...)',
    label="Max Benzonase [mmol/L]",
    title="4d: Max Benzonase per model",
    ylabel="Benzonase [mmol/L]",
    color=[:steelblue :darkorange :darkgreen],
    legend=false, size=(500,380))
savefig(fig4d, "h4d_samenvatting_drie_modellen.png")
println("  Figuur opgeslagen: h4d_samenvatting_drie_modellen.png")

println("\n=== Hoofdstuk 4 voltooid ===")
