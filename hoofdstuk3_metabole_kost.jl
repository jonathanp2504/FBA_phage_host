# ============================================================
#  HOOFDSTUK 3: Metabole kost van Benzonase
#  Model 2 (groei-verlies) vs Model 3 (Benzonase in FBA)
#
#  Analyses:
#   3a. Groeisnelheid lysogeen: Model 2 vs Model 3 vs f_prod
#   3b. Benzonase tijdsreeks: beide modellen bij zelfde parameters
#   3c. Optimizer vergelijking: optimale procesparameters per model
#   3d. Gevoeligheidsanalyse f_prod
# ============================================================

using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels, OrdinaryDiffEqCore

if !isdefined(OrdinaryDiffEqCore, :DEVerbosity)
    Core.eval(OrdinaryDiffEqCore, :(const DEVerbosity = () -> true))
end

# Laad Model 3 (Benzonase in FBA) — dit is je huidige meest complexe model
include("./Bin_model_3/parameters.jl")
include("./Bin_model_3/FBA.jl")
include("./Bin_model_3/dFBA.jl")
include("optimization_weighted_B.jl")

model_path    = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn     = 2.0;  beta_deg = 0.5
K_s           = [0.0278, 0.0146, 0.0543, 0.0833]
tau           = 1.0;  b = 170.0
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

naiveModel   = loadFBAmodel(model_path)
lysogenModel = addBenzonase!(loadFBAmodel(model_path), benz_stoich)
all_ex_ids   = [id for id in keys(naiveModel.reactions) if startswith(id, "R_EX_")]
naiveFba     = buildFbaCache(naiveModel,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba   = buildFbaCache(lysogenModel, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                              benz_id="R_BENZ_prod")

function make_p3(N0, t_inf, moi, f_prod_val=0.0015)
    Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba, lysogenFba,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        1e-10, 10.0, 5.0, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.05, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod_val)
end

# ============================================================
#  3a. Groeisnelheid lysogeen over de tijd: effect van f_prod
#  Hoe meer f_prod, hoe meer metabole last, hoe lager mu_l
# ============================================================
println("=== 3a: Groeisnelheid lysogeen vs f_prod ===")

f_prod_values = [0.0005, 0.001, 0.0015, 0.003, 0.005, 0.01]
ref_sol = run(make_p3(1e9, 2.0, 2.0, 0.0015))

# Sla mu_l op over de tijd voor elk f_prod
# mu_l is opgeslagen in p.mu_l maar verandert elke minuut
# We benaderen via de groei van de lysogene populatie
mu_l_gemiddeld = Float64[]
benz_eindwaarde = Float64[]

for fp in f_prod_values
    sol_fp = run(make_p3(1e9, 2.0, 2.0, fp))
    if SciMLBase.successful_retcode(sol_fp)
        # Schat gemiddelde groeisnelheid uit lysogene populatiegroei
        Lys_ts = [sol_fp.u[i][lind] for i in eachindex(sol_fp.u)]
        B_ts   = [sol_fp.u[i][Benzind] for i in eachindex(sol_fp.u)]
        # Bereken groeisnelheid als ln(L(t2)/L(t1))/(t2-t1) in groeifase
        idx_start = findfirst(x -> x > 1e4, Lys_ts)
        idx_mid   = isnothing(idx_start) ? 2 : min(idx_start+30, length(Lys_ts))
        if !isnothing(idx_start) && Lys_ts[idx_mid] > Lys_ts[idx_start]
            mu_est = log(Lys_ts[idx_mid]/Lys_ts[idx_start]) /
                     (sol_fp.t[idx_mid] - sol_fp.t[idx_start])
            push!(mu_l_gemiddeld, max(0.0, mu_est))
        else
            push!(mu_l_gemiddeld, 0.0)
        end
        push!(benz_eindwaarde, maximum(B_ts))
        println("  f_prod=$fp: benz=$(round(benz_eindwaarde[end],sigdigits=2)) mmol/L")
    end
end

pa = plot(f_prod_values, benz_eindwaarde,
    label="Max Benzonase", marker=:circle, lw=2, color=:green,
    xlabel="f_prod [-]", ylabel="Max Benzonase [mmol/L]",
    title="3a: Benzonase productie vs metabole last")
pb = plot(f_prod_values, mu_l_gemiddeld,
    label="Geschatte µ_l", marker=:square, lw=2, color=:orange,
    xlabel="f_prod [-]", ylabel="µ_l [h⁻¹]",
    title="3a: Groeisnelheid lysogeen vs f_prod")
fig3a = plot(pa, pb, layout=(1,2), size=(900,380))
savefig(fig3a, "h3a_fprod_effect.png")
println("  Figuur opgeslagen: h3a_fprod_effect.png")

# ============================================================
#  3b. Benzonase tijdsreeks: Model 2 (groei-verlies) vs Model 3 (FBA)
#  Zelfde f_prod, zelfde beginparameters
# ============================================================
println("\n=== 3b: Benzonase tijdsreeks Model 2 vs Model 3 ===")

# Model 3 simulatie (Benzonase in FBA)
sol_m3 = run(make_p3(1e9, 2.0, 2.0, 0.0015))
t_m3   = sol_m3.t
B_m3   = [sol_m3.u[i][Benzind] for i in eachindex(sol_m3.u)]
Lys_m3 = [sol_m3.u[i][lind]   for i in eachindex(sol_m3.u)]
N_m3   = [sol_m3.u[i][Nind]   for i in eachindex(sol_m3.u)]

# Model 2 simulatie (groei-verlies, laad aparte parameters)
# AANPASSING NODIG: laad hier Model 2 parameters en run
# Als placeholder gebruiken we Model 3 met Y_benz berekend
# zodat het vergelijkbaar is — pas aan naar jouw Model 2 bestanden
println("  Model 3 max Benzonase: $(round(maximum(B_m3), sigdigits=3)) mmol/L")
println("  (Laad hier Model 2 voor directe vergelijking)")

fig3b = plot(t_m3, B_m3,
    label="Model 3 (Benzonase in FBA)", color=:green, lw=2,
    xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    title="3b: Benzonase productie tijdsreeks",
    legend=:topleft)
# Voeg Model 2 lijn toe zodra je de Model 2 simulatie hebt:
# plot!(fig3b, t_m2, B_m2, label="Model 2 (groei-verlies)",
#       color=:blue, lw=2, linestyle=:dash)

pc = plot(t_m3, [N_m3 Lys_m3],
    label=["Naïef N" "Lysogeen l"],
    color=[:blue :green], lw=2,
    yscale=:log10, ylims=(1,:auto),
    title="3b: Populatiedynamica Model 3",
    xlabel="t [h]", ylabel="cellen/L")
fig3b_full = plot(fig3b, pc, layout=(1,2), size=(900,380))
savefig(fig3b_full, "h3b_benzonase_tijdsreeks.png")
println("  Figuur opgeslagen: h3b_benzonase_tijdsreeks.png")

# ============================================================
#  3c. Optimizer vergelijking: optimale procesparameters Model 3
# ============================================================
println("\n=== 3c: Optimalisatie Model 3 ===")

p_base = make_p3(1e9, 2.0, 2.0)
opt    = run_optimization_B(p_base)

println("  Optimaal t_inf      : $(round(opt.best_t_inf,   digits=2)) h")
println("  Optimale biomassa   : $(round(opt.best_biomass, sigdigits=2)) cellen/L")
println("  Gewogen score       : $(round(opt.best_score,   digits=5))")

p_opt   = make_p3(opt.best_biomass, opt.best_t_inf, FIXED_MOI_B)
sol_opt = run(p_opt)
B_opt   = [sol_opt.u[i][Benzind] for i in eachindex(sol_opt.u)]
Pf_opt  = [sol_opt.u[i][Pfind]  for i in eachindex(sol_opt.u)]
N_opt   = [sol_opt.u[i][Nind]   for i in eachindex(sol_opt.u)]
PN_opt  = [N_opt[i] > 1.0 ? Pf_opt[i]/N_opt[i] : 0.0 for i in eachindex(sol_opt.t)]

pd = plot(sol_opt.t, B_opt,
    label="Benzonase (optimaal)", color=:green, lw=2, fill=(0,0.2,:green),
    xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    title="3c: Optimale Benzonase productie (Model 3)")
pe = plot(sol_opt.t, PN_opt,
    label="P/N verhouding", color=:red, lw=2,
    yscale=:log10, xlabel="t [h]", ylabel="P/N [-]",
    title="3c: P/N verhouding over de tijd (Model 3)")
fig3c = plot(pd, pe, layout=(1,2), size=(900,380))
savefig(fig3c, "h3c_optimalisatie_model3.png")
println("  Figuur opgeslagen: h3c_optimalisatie_model3.png")

# ============================================================
#  3d. Gevoeligheidsanalyse f_prod op optimale score
# ============================================================
println("\n=== 3d: Gevoeligheidsanalyse f_prod ===")

scores_fp = Float64[]
for fp in f_prod_values
    p_fp  = make_p3(1e9, 2.0, 2.0, fp)
    sol_fp = run(p_fp)
    if SciMLBase.successful_retcode(sol_fp)
        B_max = maximum([sol_fp.u[i][Benzind] for i in eachindex(sol_fp.u)])
        Pf_f  = sol_fp.u[end][Pfind]
        N_f   = sol_fp.u[end][Nind]
        PN_f  = N_f > 1.0 ? Pf_f/N_f : 1e12
        score = 0.8*(B_max/BENZ_REF_B) - 0.2*(PN_f/PN_REF_B)
        push!(scores_fp, score)
        println("  f_prod=$fp: score=$(round(score,digits=4))")
    else
        push!(scores_fp, -Inf)
    end
end

fig3d = plot(f_prod_values, scores_fp,
    marker=:diamond, lw=2, color=:purple,
    xlabel="f_prod [-]",
    ylabel="Gewogen score (0.8B - 0.2PN)",
    title="3d: Gewogen score als functie van f_prod",
    label="Gewogen score", legend=false)
vline!(fig3d, [0.0015], label="Standaard f_prod=0.0015",
    color=:black, linestyle=:dash)
savefig(fig3d, "h3d_fprod_gevoeligheid.png")
println("  Figuur opgeslagen: h3d_fprod_gevoeligheid.png")

println("\n=== Hoofdstuk 3 voltooid ===")
