# ============================================================
#  HOOFDSTUK 4: Burst size — vast, lineair of metabolisch?
#  Model 3 (vast b) vs Model 4 (lineair b(µ)) vs Model 5 (FBA b)
# ============================================================
include("Model3.jl")
include("Model4.jl")
include("Model5.jl")

using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels
# ============================================================
#  Gedeelde parameters
# ============================================================
model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 3.75, 0.0, 4.0]
E_coli_cellDW  = 1e-12
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 6.0e-12;  duration = 20.0
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
N0=1e9; t_inf=2.0; moi=2.0; f_prod=0.0015

# FBA caches
naiveModel3   = Model3.loadFBAmodel(model_path)
lysogenModel3 = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
all_ex_ids    = [id for id in keys(naiveModel3.reactions) if startswith(id, "R_EX_")]
naiveFba3     = Model3.buildFbaCache(naiveModel3,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3   = Model3.buildFbaCache(lysogenModel3, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                      benz_id="R_BENZ_prod")

naiveModel4   = Model4.loadFBAmodel(model_path)
lysogenModel4 = Model4.addBenzonase!(Model4.loadFBAmodel(model_path), Model4.benz_stoich)
naiveFba4     = Model4.buildFbaCache(naiveModel4,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba4   = Model4.buildFbaCache(lysogenModel4, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                      benz_id="R_BENZ_prod")

naiveModel5   = Model5.loadFBAmodel(model_path)
lysogenModel5 = Model5.addBenzonase!(Model5.loadFBAmodel(model_path), Model5.benz_stoich)
lyticModel5   = Model5.addPhage!(Model5.loadFBAmodel(model_path), Model5.phage_stoich)
naiveFba5     = Model5.buildFbaCache(naiveModel5,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba5   = Model5.buildFbaCache(lysogenModel5, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                      benz_id="R_BENZ_prod")
lyticFba5     = Model5.buildLyticFbaCache(lyticModel5, exchange_ids,
                                           "R_BIOMASS_Ec_iJO1366_core_53p95M", "R_PHAGE_prod")

# ============================================================
#  Simuleer Model 3 (vast b=170)
# ============================================================
println("=== Simulatie Model 3 (vast b=170) ===")

p3 = Model3.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, 170.0, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba3, lysogenFba3,
    0.0, zeros(4), 0.0, zeros(4), 0.0,
    1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
    "R_BENZ_prod", 0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)
sol3 = Model3.run(p3)

t3   = sol3.t
B3   = [sol3.u[i][Model3.Benzind] for i in eachindex(sol3.u)]
Pf3  = [max(sol3.u[i][Model3.Pfind], 1.0)  for i in eachindex(sol3.u)]
Lys3 = [sol3.u[i][Model3.lind]   for i in eachindex(sol3.u)]
N3   = [sol3.u[i][Model3.Nind]   for i in eachindex(sol3.u)]
println("  Max Benzonase: $(round(maximum(B3), sigdigits=3)) mmol/L")

# ============================================================
#  Simuleer Model 4 (lineaire b(µ) = b0 + kb*µN)
#  Sla mu_N op via een wrapper om de burst size tijdsreeks te reconstrueren
# ============================================================
println("=== Simulatie Model 4 (lineaire burst size) ===")

b0 = 50.0; kb = 95.2
p4 = Model4.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, b0, kb, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba4, lysogenFba4,
    0.0, zeros(4), 0.0, zeros(4), 0.0,
    1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
    "R_BENZ_prod", 0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)

# Sla mu_N op tijdens simulatie voor burst size reconstructie
mu_N_tijds_m4 = Float64[]
mu_N_t_m4     = Float64[]

original_fba4 = Model4.fbaUpdate!
function fbaUpdate_track!(u, p)
    Model4.fbaUpdate!(u, p)
    push!(mu_N_tijds_m4, p.mu_N)
    push!(mu_N_t_m4, NaN)  # tijdstip niet direct beschikbaar hier
end

sol4 = Model4.run(p4)

t4   = sol4.t
B4   = [sol4.u[i][Model4.Benzind] for i in eachindex(sol4.u)]
Pf4  = [max(sol4.u[i][Model4.Pfind], 1.0)  for i in eachindex(sol4.u)]
Lys4 = [sol4.u[i][Model4.lind]   for i in eachindex(sol4.u)]
N4   = [sol4.u[i][Model4.Nind]   for i in eachindex(sol4.u)]

# Benadering burst size: gebruik gemiddelde mu_N uit de simulatie
# De MOI toestandsvariabele loopt van 0 naar 1, en mu_N is gebonden aan de groei
# We gebruiken de groei van N4 om mu_N te reconstrueren
# b(t) = b0 + kb * mu_N(t)
# mu_N(t) ≈ d(ln N)/dt -> numerieke afgeleid
mu_N_recon = zeros(length(t4))
for i in 2:length(t4)-1
    dt = t4[i+1] - t4[i-1]
    if N4[i] > 1.0 && N4[i-1] > 1.0 && dt > 0
        mu_N_recon[i] = max(0.0, log(max(N4[i+1],1.0)/max(N4[i-1],1.0)) / dt)
    end
end
mu_N_recon[1] = mu_N_recon[2]
mu_N_recon[end] = mu_N_recon[end-1]
b_tijds_m4 = b0 .+ kb .* mu_N_recon

println("  Max Benzonase: $(round(maximum(B4), sigdigits=3)) mmol/L")

# ============================================================
#  Simuleer Model 5 (FBA burst size)
# ============================================================
println("=== Simulatie Model 5 (FBA burst size) ===")

p5 = Model5.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
    tau, E_coli_cellDW, MW_values, h_release,
    "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
    naiveFba5, lysogenFba5, lyticFba5,
    0.0, zeros(4), 0.0, zeros(4), 0.0, 0.0,
    1e-10, 10.0, 5.0, 0.01, t_inf, moi*N0,
    "R_BENZ_prod", "R_PHAGE_prod",
    0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_prod)
sol5 = Model5.run(p5)

t5   = sol5.t
B5   = [sol5.u[i][Model5.Benzind] for i in eachindex(sol5.u)]
Pf5  = [max(sol5.u[i][Model5.Pfind], 1.0)  for i in eachindex(sol5.u)]
Lys5 = [sol5.u[i][Model5.lind]   for i in eachindex(sol5.u)]
N5   = [sol5.u[i][Model5.Nind]   for i in eachindex(sol5.u)]
println("  Max Benzonase: $(round(maximum(B5), sigdigits=3)) mmol/L")

# ============================================================
#  4a. Burst size tijdsverloop
# ============================================================
println("\n=== 4a: Burst size tijdsverloop ===")

b_tijds_m3 = fill(170.0, length(t3))

fig4a = plot(t3, b_tijds_m3, label="Model 3: vast b=170",
    color=:steelblue, lw=2, linestyle=:dash,
    xlabel="t [h]", ylabel="Burst size b [fagen/cel]",
    title="4a: Burst size over de tijd", ylims=(0, 250), legend=:topright)
plot!(fig4a, t4, b_tijds_m4, label="Model 4: lineair b(µ)", color=:darkorange, lw=2)
annotate!(fig4a, duration*0.7, 180, text("Model 5: b uit FBA\n(zie tekst)", 8, :darkgreen))
savefig(fig4a, "h4a_burst_size_tijdsverloop.png")
println("  Figuur opgeslagen: h4a_burst_size_tijdsverloop.png")

# ============================================================
#  4b. Benzonase productie
# ============================================================
println("\n=== 4b: Benzonase vergelijking ===")

fig4b = plot(t3, B3, label="Model 3 (vast b)", color=:steelblue, lw=2)
plot!(fig4b, t4, B4, label="Model 4 (lineair b)", color=:darkorange, lw=2, linestyle=:dash)
plot!(fig4b, t5, B5, label="Model 5 (FBA b)", color=:darkgreen, lw=2, linestyle=:dot)
xlabel!(fig4b, "t [h]"); ylabel!(fig4b, "Benzonase [mmol/L]")
title!(fig4b, "4b: Benzonase productie — Model 3 vs 4 vs 5")
savefig(fig4b, "h4b_benzonase_drie_modellen.png")
println("  Figuur opgeslagen: h4b_benzonase_drie_modellen.png")

# ============================================================
#  4c. Faagpopulatiedynamica
# ============================================================
println("\n=== 4c: Faagpopulatiedynamica ===")

fig4c = plot(t3, Pf3, label="Model 3 (vast b)", color=:steelblue, lw=2,
    yscale=:log10, ylims=(1,:auto), xlabel="t [h]", ylabel="Vrije fagen/L",
    title="4c: Faagpopulatie — drie modellen")
plot!(fig4c, t4, Pf4, label="Model 4 (lineair b)", color=:darkorange, lw=2, linestyle=:dash)
plot!(fig4c, t5, Pf5, label="Model 5 (FBA b)", color=:darkgreen, lw=2, linestyle=:dot)
savefig(fig4c, "h4c_fagen_drie_modellen.png")
println("  Figuur opgeslagen: h4c_fagen_drie_modellen.png")

# ============================================================
#  4d. Samenvattingstabel
# ============================================================
println("\n=== 4d: Samenvattingstabel ===")

function pn_eind(sol, Pfind_, Nind_)
    P_f = sol.u[end][Pfind_]; N_f = sol.u[end][Nind_]
    return N_f > 1.0 ? P_f/N_f : Inf
end

benz_ref = max(maximum(B3), maximum(B4), maximum(B5))
pn3 = pn_eind(sol3, Model3.Pfind, Model3.Nind)
pn4 = pn_eind(sol4, Model4.Pfind, Model4.Nind)
pn5 = pn_eind(sol5, Model5.Pfind, Model5.Nind)
pn_vals = filter(!isinf, [pn3, pn4, pn5])
pn_ref  = isempty(pn_vals) ? 1.0 : max(pn_vals...)

results = [
    ("Model 3", maximum(B3), pn3),
    ("Model 4", maximum(B4), pn4),
    ("Model 5", maximum(B5), pn5),
]

println("\n  Model   | Max Benz [mmol/L] | P/N eind  | Score (0.8B-0.2PN)")
println("  --------|-------------------|-----------|--------------------")
for (naam, benz, pn) in results
    score = 0.8*(benz/benz_ref) - 0.2*(isinf(pn) ? 1.0 : pn/pn_ref)
    @printf("  %-7s | %-17.5f | %-9.2e | %.4f\n", naam, benz, pn, score)
end

fig4d = bar(["Model 3" "Model 4" "Model 5"],
    [maximum(B3) maximum(B4) maximum(B5)],
    title="4d: Max Benzonase per model",
    ylabel="Benzonase [mmol/L]",
    color=[:steelblue :darkorange :darkgreen],
    legend=false, size=(500,380))
savefig(fig4d, "h4d_samenvatting_drie_modellen.png")
println("  Figuur opgeslagen: h4d_samenvatting_drie_modellen.png")

println("\n=== Hoofdstuk 4 voltooid ===")