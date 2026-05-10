# ============================================================
#  HOOFDSTUK 3: Metabole kost van Benzonase
#  Model 2 (groei-verlies) vs Model 3 (Benzonase in FBA)
#
#  Analyses:
#   3a. Benzonase productie en groeisnelheid vs f_prod (Model 3)
#   3b. Benzonase tijdsreeks: Model 2 vs Model 3
#   3c. Optimale run Model 3 (waarden uit run_optimalisatie_3.jl)
#   3d. Gevoeligheidsanalyse f_prod (Model 3)
#
#  Optimale waarden Model 3 (uit run_optimalisatie_3.jl):
#   Optimizer A: MOI=1.2289 | t_inf=0.5h  | biomassa=1e9   | score=1.319748
#   Optimizer B: MOI=2.0    | t_inf=1.75h | biomassa=1.21e8 | score=1.30051
# ============================================================
include("Model2.jl")
include("Model3.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels, OrdinaryDiffEqCore


# Referentiewaarden voor normalisatie (uit run_kalibratie_3.jl)
const BENZ_REF_3 = 0.001
const PN_REF_3   = 2.0e10

# ============================================================
#  Gedeelde parameters
# ============================================================
model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0;  b = 170.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 3.75, 0.0, 4.0]
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

# FBA caches Model 2
naiveModel2    = Model2.loadFBAmodel(model_path)
all_ex_ids     = [id for id in keys(naiveModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2      = Model2.buildFbaCache(naiveModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2    = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                       "R_BIOMASS_Ec_iJO1366_core_53p95M")

# FBA caches Model 3
naiveModel3    = Model3.loadFBAmodel(model_path)
lysogenModel3  = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
naiveFba3      = Model3.buildFbaCache(naiveModel3,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3    = Model3.buildFbaCache(lysogenModel3, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                       benz_id="R_BENZ_prod")

function make_p2(N0, t_inf, moi, f_prod_val=0.0015)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        1e-10, 10.0, 5.0, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0,
        mu_max_vector, e_max_vector, f_prod_val, 0.05)
end

function make_p3(N0, t_inf, moi, f_prod_val=0.0015)
    Model3.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        1e-10, 10.0, 5.0, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.05, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod_val)
end

# ============================================================
#  3a. Benzonase productie en groeisnelheid vs f_prod (Model 3)
# ============================================================
println("=== 3a: Groeisnelheid lysogeen vs f_prod ===")

f_prod_values   = [0.0005, 0.001, 0.0015, 0.003, 0.005, 0.01]
mu_l_gemiddeld  = Float64[]
benz_eindwaarde = Float64[]

for fp in f_prod_values
    sol_fp = Model3.run(make_p3(1e9, 2.0, 2.0, fp))
    if SciMLBase.successful_retcode(sol_fp)
        Lys_ts = [sol_fp.u[i][Model3.lind]   for i in eachindex(sol_fp.u)]
        B_ts   = [sol_fp.u[i][Model3.Benzind] for i in eachindex(sol_fp.u)]
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

pa = plot(f_prod_values, benz_eindwaarde, label="Max Benzonase (Model 3)",
    marker=:circle, lw=2, color=:green,
    xlabel="f_prod [-]", ylabel="Max Benzonase [mmol/L]",
    title="3a: Benzonase productie vs metabole last")
pb = plot(f_prod_values, mu_l_gemiddeld, label="Geschatte µ_l (Model 3)",
    marker=:square, lw=2, color=:orange,
    xlabel="f_prod [-]", ylabel="µ_l [h⁻¹]",
    title="3a: Groeisnelheid lysogeen vs f_prod")
fig3a = plot(pa, pb, layout=(1,2), size=(900,380))
savefig(fig3a, "h3a_fprod_effect.png")
println("  Figuur opgeslagen: h3a_fprod_effect.png")

# ============================================================
#  3b. Benzonase tijdsreeks: Model 2 vs Model 3
#  Zelfde parameters zodat de vergelijking eerlijk is
# ============================================================
println("\n=== 3b: Benzonase tijdsreeks Model 2 vs Model 3 ===")

sol_m2 = Model2.run(make_p2(1e9, 2.0, 2.0, 0.0015))
sol_m3 = Model3.run(make_p3(1e9, 2.0, 2.0, 0.0015))

t_m2   = sol_m2.t
B_m2   = [sol_m2.u[i][Model2.Benzind] for i in eachindex(sol_m2.u)]
t_m3   = sol_m3.t
B_m3   = [sol_m3.u[i][Model3.Benzind] for i in eachindex(sol_m3.u)]
Lys_m3 = [sol_m3.u[i][Model3.lind]   for i in eachindex(sol_m3.u)]
N_m3   = [sol_m3.u[i][Model3.Nind]   for i in eachindex(sol_m3.u)]

println("  Model 2 max Benzonase: $(round(maximum(B_m2), sigdigits=3)) mmol/L")
println("  Model 3 max Benzonase: $(round(maximum(B_m3), sigdigits=3)) mmol/L")

fig3b = plot(t_m2, B_m2, label="Model 2 (groei-verlies)",
    color=:steelblue, lw=2, linestyle=:dash)
plot!(fig3b, t_m3, B_m3, label="Model 3 (Benzonase in FBA)",
    color=:green, lw=2,
    xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    title="3b: Benzonase productie tijdsreeks", legend=:topleft)

pc = plot(t_m3, [N_m3 Lys_m3], label=["Naïef N" "Lysogeen l"],
    color=[:blue :green], lw=2, yscale=:log10, ylims=(1,:auto),
    title="3b: Populatiedynamica Model 3", xlabel="t [h]", ylabel="cellen/L")
fig3b_full = plot(fig3b, pc, layout=(1,2), size=(900,380))
savefig(fig3b_full, "h3b_benzonase_tijdsreeks.png")
println("  Figuur opgeslagen: h3b_benzonase_tijdsreeks.png")

# ============================================================
#  3c. Optimale run Model 3
#  Waarden uit run_optimalisatie_3.jl — geen optimizer aanroep
#
#  Optimizer A (beste MOI + t_inf):
#    MOI=1.2289 | t_inf=0.5h | biomassa=1e9 | score=1.319748
#  Optimizer B (beste t_inf + biomassa, MOI=2.0):
#    MOI=2.0 | t_inf=1.75h | biomassa=1.21e8 | score=1.30051
# ============================================================
println("\n=== 3c: Optimale run Model 3 ===")

# Optimizer A resultaat
opt_A_moi      = 1.2289
opt_A_t_inf    = 0.5
opt_A_biomassa = 1e9
opt_A_score    = 1.319748

# Optimizer B resultaat
opt_B_moi      = 2.0
opt_B_t_inf    = 1.75
opt_B_biomassa = 1.21e8
opt_B_score    = 1.30051

println("  Optimizer A: MOI=$opt_A_moi | t_inf=$opt_A_t_inf h | biomassa=$opt_A_biomassa | score=$opt_A_score")
println("  Optimizer B: MOI=$opt_B_moi | t_inf=$opt_B_t_inf h | biomassa=$opt_B_biomassa | score=$opt_B_score")
println("  → Beste resultaat: Optimizer A (hogere score)")

# Simuleer optimale run (Optimizer A heeft hoogste score)
p_opt_A   = make_p3(opt_A_biomassa, opt_A_t_inf, opt_A_moi)
sol_opt_A = Model3.run(p_opt_A)
B_opt_A   = [sol_opt_A.u[i][Model3.Benzind] for i in eachindex(sol_opt_A.u)]
Pf_opt_A  = [sol_opt_A.u[i][Model3.Pfind]  for i in eachindex(sol_opt_A.u)]
N_opt_A   = [sol_opt_A.u[i][Model3.Nind]   for i in eachindex(sol_opt_A.u)]
PN_opt_A  = [N_opt_A[i] > 1.0 ? Pf_opt_A[i]/N_opt_A[i] : 0.0
             for i in eachindex(sol_opt_A.t)]

# Simuleer ook Optimizer B voor vergelijking
p_opt_B   = make_p3(opt_B_biomassa, opt_B_t_inf, opt_B_moi)
sol_opt_B = Model3.run(p_opt_B)
B_opt_B   = [sol_opt_B.u[i][Model3.Benzind] for i in eachindex(sol_opt_B.u)]

println("  Max Benzonase (Opt A): $(round(maximum(B_opt_A), sigdigits=3)) mmol/L")
println("  Max Benzonase (Opt B): $(round(maximum(B_opt_B), sigdigits=3)) mmol/L")

pd = plot(sol_opt_A.t, B_opt_A, label="Optimizer A (MOI=1.23, t_inf=0.5h)",
    color=:green, lw=2, fill=(0,0.2,:green))
plot!(pd, sol_opt_B.t, B_opt_B, label="Optimizer B (MOI=2.0, t_inf=1.75h)",
    color=:darkgreen, lw=2, linestyle=:dash,
    xlabel="t [h]", ylabel="Benzonase [mmol/L]",
    title="3c: Optimale Benzonase productie (Model 3)")

pe = plot(sol_opt_A.t, PN_opt_A, label="P/N (Optimizer A)", color=:red, lw=2,
    yscale=:log10, xlabel="t [h]", ylabel="P/N [-]",
    title="3c: P/N verhouding (Model 3, Optimizer A)")
fig3c = plot(pd, pe, layout=(1,2), size=(900,380))
savefig(fig3c, "h3c_optimalisatie_model3.png")
println("  Figuur opgeslagen: h3c_optimalisatie_model3.png")

# ============================================================
#  3d. Gevoeligheidsanalyse f_prod (Model 3)
# ============================================================
println("\n=== 3d: Gevoeligheidsanalyse f_prod ===")

scores_fp = Float64[]
for fp in f_prod_values
    sol_fp = Model3.run(make_p3(1e9, 2.0, 2.0, fp))
    if SciMLBase.successful_retcode(sol_fp)
        B_max = maximum([sol_fp.u[i][Model3.Benzind] for i in eachindex(sol_fp.u)])
        Pf_f  = sol_fp.u[end][Model3.Pfind]
        N_f   = sol_fp.u[end][Model3.Nind]
        PN_f  = N_f > 1.0 ? Pf_f/N_f : 1e12
        score = 0.8*(B_max/BENZ_REF_3) - 0.05*(PN_f/PN_REF_3)
        push!(scores_fp, score)
        println("  f_prod=$fp: score=$(round(score,digits=4))")
    else
        push!(scores_fp, -Inf)
    end
end

fig3d = plot(f_prod_values, scores_fp, marker=:diamond, lw=2, color=:purple,
    xlabel="f_prod [-]", ylabel="Gewogen score (0.8B - 0.05PN)",
    title="3d: Gewogen score als functie van f_prod (Model 3)", legend=false)
vline!(fig3d, [0.0015], label="Standaard f_prod=0.0015",
    color=:black, linestyle=:dash)
savefig(fig3d, "h3d_fprod_gevoeligheid.png")
println("  Figuur opgeslagen: h3d_fprod_gevoeligheid.png")

println("\n=== Hoofdstuk 3 voltooid ===")