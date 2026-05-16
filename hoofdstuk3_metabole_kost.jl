# ============================================================
#  HOOFDSTUK 3: Metabole kost van Benzonase
#  Model 2 (groei-verlies) vs Model 3 (Benzonase in FBA)
#
#  Analyses:
#   3a. Werkelijke groeisnelheid lysogeen vs f_prod (Model 3)
#   3b. Benzonase tijdsreeks: Model 2 vs Model 3
#   3c. Optimale run Model 3
#   3d. Gevoeligheidsanalyse f_prod (Model 3)
#   3e. MOI-effect op Benzonase: Model 2 en Model 3 naast elkaar
#   3f. Sensitivity analysis alpha_syn
# ============================================================
include("Model2.jl")
include("Model3.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

const BENZ_REF_3 = 0.001
const PN_REF_3   = 2.0e10

model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 1.32;  b = 170.0
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [0.0, 2.26, 0.0, 10.0]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 1.71e-12;  duration = 20.0
mu_max_vector  = [0.76, 0.76, 1.10, 0.30]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW  = 1.0e-12

naiveModel2   = Model2.loadFBAmodel(model_path)
all_ex_ids    = [id for id in keys(naiveModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2     = Model2.buildFbaCache(naiveModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2   = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                      "R_BIOMASS_Ec_iJO1366_core_53p95M")
naiveModel3   = Model3.loadFBAmodel(model_path)
lysogenModel3 = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
naiveFba3     = Model3.buildFbaCache(naiveModel3,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3   = Model3.buildFbaCache(lysogenModel3, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M";
                                      benz_id="R_BENZ_prod")

function make_p2(N0, t_inf, moi, f_prod_val=0.0015)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
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
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.05, 0.1, 0.0,
        mu_max_vector, e_max_vector, f_prod_val)
end

# ============================================================
#  3a. Werkelijke groeisnelheid lysogeen vs f_prod (Model 3)
# ============================================================
println("=== 3a: Werkelijke groeisnelheid lysogeen vs f_prod ===")

k_tox_val      = 0.05
f_prod_values  = [0.0005, 0.001, 0.0015, 0.003, 0.005, 0.01]
mu_l_gemiddeld = Float64[]
benz_eindwaarde = Float64[]

for fp in f_prod_values
    sol_fp = Model3.run(make_p3(1e9, 2.0, 2.0, fp))
    if SciMLBase.successful_retcode(sol_fp)
        lind_ts  = [sol_fp.u[i][Model3.lind]   for i in eachindex(sol_fp.u)]
        B_ts     = [sol_fp.u[i][Model3.Benzind] for i in eachindex(sol_fp.u)]
        dt_vec   = diff(sol_fp.t)
        dlind_dt = diff(lind_ts) ./ dt_vec
        mu_l_ts  = [(dlind_dt[i] / max(lind_ts[i], 1.0) + k_tox_val) / (1.0 - fp)
                    for i in eachindex(dlind_dt)]
        stabiel_idx = findall(i -> lind_ts[i] > 1e8 && mu_l_ts[i] > 0.0 && mu_l_ts[i] < 2.0,
                              eachindex(mu_l_ts))
        mu_gem = isempty(stabiel_idx) ? 0.0 : mean(mu_l_ts[stabiel_idx])
        push!(mu_l_gemiddeld,  mu_gem)
        push!(benz_eindwaarde, maximum(B_ts))
        println("  f_prod=$fp: benz=$(round(benz_eindwaarde[end],sigdigits=2)) | µ_l=$(round(mu_gem,digits=4)) h⁻¹")
    end
end

pa = plot(f_prod_values, benz_eindwaarde,
    label="Max Benzonase (Model 3)",
    marker=:circle, lw=2, color=:green,
    xlabel="f_prod [-]", ylabel="Max Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
pb = plot(f_prod_values, mu_l_gemiddeld,
    label="Geschatte µ_l (Model 3)",
    marker=:square, lw=2, color=:orange,
    xlabel="f_prod [-]", ylabel="µ_l [h⁻¹]",
    legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
fig3a = plot(pa, pb, layout=(1,2), size=(900,420), margin=6Plots.mm)
savefig(fig3a, "h3a_fprod_effect.png")
println("  Figuur opgeslagen: h3a_fprod_effect.png")

# ============================================================
#  3b. Benzonase tijdsreeks: Model 2 vs Model 3
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

p3b_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p3b_benz, t_m2, B_m2, label="Model 2 (groei-verlies)",
    color=:steelblue, lw=2, linestyle=:dash)
plot!(p3b_benz, t_m3, B_m3, label="Model 3 (FBA)",
    color=:green, lw=2)

p3b_pop = plot(xlabel="t [h]", ylabel="Cellen L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p3b_pop, t_m3, [N_m3 Lys_m3],
    label=["Naïef N" "Lysogeen l"],
    color=[:blue :green], lw=2)

fig3b_full = plot(p3b_benz, p3b_pop, layout=(1,2), size=(900,430), margin=6Plots.mm)
savefig(fig3b_full, "h3b_benzonase_tijdsreeks.png")
println("  Figuur opgeslagen: h3b_benzonase_tijdsreeks.png")

# ============================================================
#  3c. Optimale run Model 3
# ============================================================
println("\n=== 3c: Optimale run Model 3 ===")

opt_A_moi=1.2289; opt_A_t_inf=0.5;  opt_A_biomassa=1e9;   opt_A_score=1.319748
opt_B_moi=2.0;    opt_B_t_inf=1.75; opt_B_biomassa=1.21e8; opt_B_score=1.30051

p_opt_A   = make_p3(opt_A_biomassa, opt_A_t_inf, opt_A_moi)
sol_opt_A = Model3.run(p_opt_A)
B_opt_A   = [sol_opt_A.u[i][Model3.Benzind] for i in eachindex(sol_opt_A.u)]
Pf_opt_A  = [sol_opt_A.u[i][Model3.Pfind]  for i in eachindex(sol_opt_A.u)]
N_opt_A   = [sol_opt_A.u[i][Model3.Nind]   for i in eachindex(sol_opt_A.u)]
PN_opt_A  = [N_opt_A[i]>1.0 ? Pf_opt_A[i]/N_opt_A[i] : 0.0
             for i in eachindex(sol_opt_A.t)]
p_opt_B   = make_p3(opt_B_biomassa, opt_B_t_inf, opt_B_moi)
sol_opt_B = Model3.run(p_opt_B)
B_opt_B   = [sol_opt_B.u[i][Model3.Benzind] for i in eachindex(sol_opt_B.u)]

pd = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(pd, sol_opt_A.t, B_opt_A,
    label="Optimizer A (MOI=1.23, t_inf=0.5h)",
    color=:green, lw=2, fill=(0,0.2,:green))
plot!(pd, sol_opt_B.t, B_opt_B,
    label="Optimizer B (MOI=2.0, t_inf=1.75h)",
    color=:darkgreen, lw=2, linestyle=:dash)

pe = plot(xlabel="t [h]", ylabel="P/N [-]",
    yscale=:log10, legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(pe, sol_opt_A.t, max.(PN_opt_A, 1e-10),
    label="P/N (Optimizer A)", color=:red, lw=2)

fig3c = plot(pd, pe, layout=(1,2), size=(900,430), margin=6Plots.mm)
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
    else
        push!(scores_fp, -Inf)
    end
end

fig3d = plot(f_prod_values, scores_fp,
    marker=:diamond, lw=2, color=:purple,
    xlabel="f_prod [-]", ylabel="Gewogen score (0.8B − 0.05PN)",
    legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
vline!(fig3d, [0.0015], color=:black, lw=1.5, linestyle=:dash,
    label="f_prod = 0.0015")
savefig(fig3d, "h3d_fprod_gevoeligheid.png")
println("  Figuur opgeslagen: h3d_fprod_gevoeligheid.png")

# ============================================================
#  3e. MOI-effect op Benzonase productie: Model 2 én Model 3
# ============================================================
println("\n=== 3e: MOI-effect op Benzonase productie (Model 2 en Model 3) ===")

moi_sweep_benz = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
benz_m2_per_moi  = Float64[]
benz_m3_per_moi  = Float64[]
lys_m2_per_moi   = Float64[]
lys_m3_per_moi   = Float64[]

for moi in moi_sweep_benz
    # Model 2
    sol_m2_moi = Model2.run(make_p2(1e9, 2.0, moi, 0.0015))
    if SciMLBase.successful_retcode(sol_m2_moi)
        B_m2_moi   = maximum([sol_m2_moi.u[i][Model2.Benzind] for i in eachindex(sol_m2_moi.u)])
        Lys_m2_moi = maximum([sol_m2_moi.u[i][Model2.lind]    for i in eachindex(sol_m2_moi.u)])
        push!(benz_m2_per_moi, B_m2_moi)
        push!(lys_m2_per_moi,  Lys_m2_moi)
    else
        push!(benz_m2_per_moi, 0.0); push!(lys_m2_per_moi, 0.0)
    end
    # Model 3
    sol_m3_moi = Model3.run(make_p3(1e9, 2.0, moi, 0.0015))
    if SciMLBase.successful_retcode(sol_m3_moi)
        B_m3_moi   = maximum([sol_m3_moi.u[i][Model3.Benzind] for i in eachindex(sol_m3_moi.u)])
        Lys_m3_moi = maximum([sol_m3_moi.u[i][Model3.lind]    for i in eachindex(sol_m3_moi.u)])
        push!(benz_m3_per_moi, B_m3_moi)
        push!(lys_m3_per_moi,  Lys_m3_moi)
    else
        push!(benz_m3_per_moi, 0.0); push!(lys_m3_per_moi, 0.0)
    end
    println("  MOI=$moi | M2 benz=$(round(benz_m2_per_moi[end],sigdigits=2)) | M3 benz=$(round(benz_m3_per_moi[end],sigdigits=2))")
end

p3e_benz = plot(xlabel="MOI [-]", ylabel="Max Benzonase [mmol L⁻¹]",
    xscale=:log10, legend=:topleft,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
plot!(p3e_benz, moi_sweep_benz, benz_m2_per_moi,
    marker=:circle, lw=2, color=:steelblue, label="Model 2 (groei-verlies)")
plot!(p3e_benz, moi_sweep_benz, benz_m3_per_moi,
    marker=:square, lw=2, color=:green,     label="Model 3 (FBA)")

p3e_lys = plot(xlabel="MOI [-]", ylabel="Max lysogene cellen L⁻¹",
    xscale=:log10, yscale=:log10, legend=:topleft,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
plot!(p3e_lys, moi_sweep_benz, max.(lys_m2_per_moi, 1.0),
    marker=:circle, lw=2, color=:steelblue, label="Model 2")
plot!(p3e_lys, moi_sweep_benz, max.(lys_m3_per_moi, 1.0),
    marker=:square, lw=2, color=:green,     label="Model 3")

fig3e = plot(p3e_benz, p3e_lys, layout=(1,2), size=(1000,430), margin=6Plots.mm)
savefig(fig3e, "h3e_moi_effect_benzonase.png")
println("  Figuur opgeslagen: h3e_moi_effect_benzonase.png")

# ============================================================
#  3f. Sensitivity analysis alpha_syn (Model 3)
#
#  alpha_syn bepaalt hoe snel enzymen voor nieuwe koolstofbronnen
#  worden aangemaakt na een substraatswitch (bv. glucose → maltose).
#  Een hogere alpha_syn verkort de diauxische vertraging en verhoogt
#  de LamB-expressie sneller na glucosedepletie.
#  Hier testen we het effect op: (1) Benzonase productie,
#  (2) lysogene populatiegrootte, (3) enzymkinetiek (e_mal verloop)
# ============================================================
println("\n=== 3f: Sensitivity analysis alpha_syn (Model 3) ===")

alpha_syn_waarden = [0.5, 1.0, 2.0, 4.0]   # standaard = 2.0
kleuren_3f        = [:darkblue :steelblue :green :darkorange]
labels_3f         = ["α_syn=0.5" "α_syn=1.0" "α_syn=2.0 (ref)" "α_syn=4.0"]

# We bouwen aparte FBA-caches voor elke alpha_syn omdat e_max afhangt van alpha_syn
function make_p3_asyn(N0, t_inf, moi, asyn, f_prod_val=0.0015)
    e_max_asyn = (asyn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
    Model3.Parameters(duration, N0, asyn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.05, 0.1, 0.0,
        e_max_asyn, e_max_asyn, f_prod_val)
end

# Tijdsreeksen voor alle alpha_syn waarden
p3f_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
p3f_lys  = plot(xlabel="t [h]", ylabel="Lysogene cellen L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
p3f_emal = plot(xlabel="t [h]", ylabel="e_mal (LamB proxy) [-]",
    legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)

benz_asyn    = Float64[]
t_max_asyn   = Float64[]

for (idx, asyn) in enumerate(alpha_syn_waarden)
    sol_as = Model3.run(make_p3_asyn(1e9, 2.0, 2.0, asyn))
    if SciMLBase.successful_retcode(sol_as)
        t_as   = sol_as.t
        B_as   = [sol_as.u[i][Model3.Benzind] for i in eachindex(sol_as.u)]
        Lys_as = max.([sol_as.u[i][Model3.lind]  for i in eachindex(sol_as.u)], 1.0)
        emal_as= [sol_as.u[i][Model3.Eind[2]]    for i in eachindex(sol_as.u)]
        i_max  = argmax(B_as)
        push!(benz_asyn,  B_as[i_max])
        push!(t_max_asyn, t_as[i_max])
        println("  α_syn=$asyn: max_benz=$(round(B_as[i_max],sigdigits=3)) | t_max=$(round(t_as[i_max],digits=2)) h")
        plot!(p3f_benz, t_as, B_as,   label=labels_3f[idx], color=kleuren_3f[idx], lw=2)
        plot!(p3f_lys,  t_as, Lys_as, label=labels_3f[idx], color=kleuren_3f[idx], lw=2)
        plot!(p3f_emal, t_as, emal_as,label=labels_3f[idx], color=kleuren_3f[idx], lw=2)
    end
end

# Samenvatting: max Benzonase en tijdstip per alpha_syn
p3f_bar_benz = bar(string.(alpha_syn_waarden), benz_asyn,
    xlabel="α_syn [h⁻¹]", ylabel="Max Benzonase [mmol L⁻¹]",
    color=kleuren_3f, legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
p3f_bar_t = bar(string.(alpha_syn_waarden), t_max_asyn,
    xlabel="α_syn [h⁻¹]", ylabel="Tijdstip max Benzonase [h]",
    color=kleuren_3f, legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)

# Lay-out: bovenste rij = tijdsreeksen, onderste rij = samenvatting
fig3f_ts  = plot(p3f_benz, p3f_lys, p3f_emal,
    layout=(1,3), size=(1300,430), margin=6Plots.mm)
fig3f_sum = plot(p3f_bar_benz, p3f_bar_t,
    layout=(1,2), size=(900,400), margin=6Plots.mm)

savefig(fig3f_ts,  "h3f_alphasyn_tijdsreeksen.png")
savefig(fig3f_sum, "h3f_alphasyn_samenvatting.png")
println("  Figuren opgeslagen: h3f_alphasyn_tijdsreeksen.png | h3f_alphasyn_samenvatting.png")

println("\n=== Hoofdstuk 3 voltooid ===")