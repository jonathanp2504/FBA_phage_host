# ============================================================
#  HOOFDSTUK 4: Burst size — fixed, linear or metabolic?
# ============================================================
include("Model3.jl"); include("Model4.jl"); include("Model5.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

model_path    = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn     = 2.0;  beta_deg = 0.5
K_s           = [0.0061, 9.4e-4, 0.0543, 8.33]
tau           = 1.32
p_pref        = [0.8925, 0.08925, 0.008925, 0.008925]
V_max         = [0.0, 2.26, 0.0, 10.0]
exchange_ids  = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids= ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                 "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                 "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                 "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                 "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
MW_values     = [180.16, 342.3, 92.09, 60.05]
h_release     = 1.71e-12;  duration = 40.0
mu_max_vector = [0.76, 0.76, 1.10, 0.30]
e_max_vector  = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
E_coli_cellDW = 1.0e-12
N0=1e9; t_inf=2.0; moi=2.0; f_prod=0.0015

naiveModel3   = Model3.loadFBAmodel(model_path)
lysogenModel3 = Model3.addBenzonase!(Model3.loadFBAmodel(model_path), Model3.benz_stoich)
all_ex_ids    = [id for id in keys(naiveModel3.reactions) if startswith(id, "R_EX_")]
naiveFba3     = Model3.buildFbaCache(naiveModel3,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba3   = Model3.buildFbaCache(lysogenModel3, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
naiveModel4   = Model4.loadFBAmodel(model_path)
lysogenModel4 = Model4.addBenzonase!(Model4.loadFBAmodel(model_path), Model4.benz_stoich)
naiveFba4     = Model4.buildFbaCache(naiveModel4,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba4   = Model4.buildFbaCache(lysogenModel4, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
naiveModel5   = Model5.loadFBAmodel(model_path)
lysogenModel5 = Model5.addBenzonase!(Model5.loadFBAmodel(model_path), Model5.benz_stoich)
lyticModel5   = Model5.addPhage!(Model5.loadFBAmodel(model_path), Model5.phage_stoich)
naiveFba5     = Model5.buildFbaCache(naiveModel5,   exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba5   = Model5.buildFbaCache(lysogenModel5, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M"; benz_id="R_BENZ_prod")
lyticFba5     = Model5.buildLyticFbaCache(lyticModel5, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M", "R_PHAGE_prod")

function make_p3(N0_, t_inf_, moi_, f_=f_prod)
    Model3.Parameters(duration, N0_, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, 170.0, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3, 0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf_, moi_*N0_,
        "R_BENZ_prod", 0.35, 0.1, 0.0, mu_max_vector, e_max_vector, f_)
end
function make_p4(N0_, t_inf_, moi_, f_=f_prod)
    Model4.Parameters(duration, N0_, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, 50.0, 228.6, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba4, lysogenFba4, 0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf_, moi_*N0_,
        "R_BENZ_prod", 0.35, 0.1, 0.0, mu_max_vector, e_max_vector, f_)
end
function make_p5(N0_, t_inf_, moi_, f_=f_prod)
    Model5.Parameters(duration, N0_, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba5, lysogenFba5, lyticFba5, 0.0, zeros(4), 0.0, zeros(4), 0.0, 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf_, moi_*N0_,
        "R_BENZ_prod", "R_PHAGE_prod", 0.35, 0.1, 0.0, mu_max_vector, e_max_vector, f_)
end

println("=== Reference simulations ===")
sol3 = Model3.run(make_p3(N0, t_inf, moi))
sol4 = Model4.run(make_p4(N0, t_inf, moi))
sol5 = Model5.run(make_p5(N0, t_inf, moi))

t3 = sol3.t; B3  = [sol3.u[i][Model3.Benzind] for i in eachindex(sol3.u)]
              Pf3 = [max(sol3.u[i][Model3.Pfind], 1.0) for i in eachindex(sol3.u)]
t4 = sol4.t; B4  = [sol4.u[i][Model4.Benzind] for i in eachindex(sol4.u)]
              Pf4 = [max(sol4.u[i][Model4.Pfind], 1.0) for i in eachindex(sol4.u)]
              N4  = [sol4.u[i][Model4.Nind] for i in eachindex(sol4.u)]
t5 = sol5.t; B5  = [sol5.u[i][Model5.Benzind] for i in eachindex(sol5.u)]
              Pf5 = [max(sol5.u[i][Model5.Pfind], 1.0) for i in eachindex(sol5.u)]

b_tijds_m3 = fill(170.0, length(t3))
mu_N_recon = zeros(length(t4))
for i in 2:length(t4)-1
    dt = t4[i+1] - t4[i-1]
    if N4[i] > 1.0 && dt > 0
        mu_N_recon[i] = max(0.0, log(max(N4[i+1],1.0)/max(N4[i-1],1.0)) / dt)
    end
end
mu_N_recon[1] = mu_N_recon[2]; mu_N_recon[end] = mu_N_recon[end-1]
b_tijds_m4 = 50.0 .+ 228.6 .* mu_N_recon
b_tijds_m5 = [sol5.u[i][Model5.Bind] for i in eachindex(sol5.u)]

# ============================================================
#  4a. Burst size until maltose depletion + maltose per model
# ============================================================
println("\n=== 4a: Burst size until maltose depletion ===")

function find_malt_depletion(sol, Sind_)
    for i in eachindex(sol.u)
        if sol.u[i][Sind_[2]] < 0.01 && sol.t[i] > 0.5
            return sol.t[i]
        end
    end
    return sol.t[end]
end

t_malt3 = find_malt_depletion(sol3, Model3.Sind)
t_malt4 = find_malt_depletion(sol4, Model4.Sind)
t_malt5 = find_malt_depletion(sol5, Model5.Sind)
t_malt_depl = minimum([t_malt3, t_malt4, t_malt5])
println("  Maltose depletion: M3=$(round(t_malt3,digits=2)) | M4=$(round(t_malt4,digits=2)) | M5=$(round(t_malt5,digits=2)) h")

idx3 = findall(t3 .<= t_malt_depl)
idx4 = findall(t4 .<= t_malt_depl)
idx5 = findall(t5 .<= t_malt_depl)

p4a_burst = plot(xlabel="t [h]", ylabel="Burst size b [phages cell⁻¹]",
    ylims=(0, max(maximum(b_tijds_m4[idx4])*1.15, 250)),
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p4a_burst, t3[idx3], b_tijds_m3[idx3], label="Model 3 (fixed b=170)",  color=:steelblue,  lw=2, linestyle=:dash)
plot!(p4a_burst, t4[idx4], b_tijds_m4[idx4], label="Model 4 (linear b(µ))", color=:darkorange, lw=2)
plot!(p4a_burst, t5[idx5], b_tijds_m5[idx5], label="Model 5 (FBA b)",        color=:darkgreen,  lw=2, linestyle=:dot)
vline!(p4a_burst, [t_inf],       color=:gray,  lw=1,   linestyle=:dash, label="t_inf")
vline!(p4a_burst, [t_malt_depl], color=:black, lw=1.5, linestyle=:dot,  label="Maltose depleted")

S_malt3 = [sol3.u[i][Model3.Sind[2]] for i in eachindex(sol3.u)]
S_malt4 = [sol4.u[i][Model4.Sind[2]] for i in eachindex(sol4.u)]
S_malt5 = [sol5.u[i][Model5.Sind[2]] for i in eachindex(sol5.u)]

p4a_malt = plot(xlabel="t [h]", ylabel="Maltose [mmol L⁻¹]",
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p4a_malt, t3, S_malt3, label="Model 3", color=:steelblue,  lw=2)
plot!(p4a_malt, t4, S_malt4, label="Model 4", color=:darkorange, lw=2, linestyle=:dash)
plot!(p4a_malt, t5, S_malt5, label="Model 5", color=:darkgreen,  lw=2, linestyle=:dot)
vline!(p4a_malt, [t_inf], color=:gray, lw=1, linestyle=:dash, label="t_inf")

fig4a = plot(p4a_burst, p4a_malt, layout=(1,2), size=(1100,430), margin=6Plots.mm)
savefig(fig4a, "h4a_burst_size_tijdsverloop.png")
println("  Saved: h4a_burst_size_tijdsverloop.png")

# ============================================================
#  4b. Benzonase comparison: 12 panels
# ============================================================
println("\n=== 4b: Benzonase comparison 12 panels ===")

moi_waarden  = [0.1, 0.5, 2.0]
tinf_waarden = [1.0, 2.0, 5.0]
N0_waarden   = [1e7, 1e8, 1e9]
kleur_m3 = :steelblue; kleur_m4 = :darkorange; kleur_m5 = :darkgreen

panels_4b = []
for moi_j in moi_waarden
    for tinf_k in tinf_waarden
        lbl = "MOI=$(moi_j), t_inf=$(tinf_k)h"
        p = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
            legend=:outertopright,
            tickfontsize=9, guidefontsize=10, legendfontsize=8,
            bottom_margin=6Plots.mm, left_margin=8Plots.mm,
            title=lbl, titlefontsize=9)

        s3 = Model3.run(make_p3(1e9, tinf_k, moi_j))
        s4 = Model4.run(make_p4(1e9, tinf_k, moi_j))
        s5 = Model5.run(make_p5(1e9, tinf_k, moi_j))

        SciMLBase.successful_retcode(s3) && plot!(p, s3.t,
            [s3.u[i][Model3.Benzind] for i in eachindex(s3.u)],
            label="M3", color=kleur_m3, lw=2)
        SciMLBase.successful_retcode(s4) && plot!(p, s4.t,
            [s4.u[i][Model4.Benzind] for i in eachindex(s4.u)],
            label="M4", color=kleur_m4, lw=2, linestyle=:dash)
        SciMLBase.successful_retcode(s5) && plot!(p, s5.t,
            [s5.u[i][Model5.Benzind] for i in eachindex(s5.u)],
            label="M5", color=kleur_m5, lw=2, linestyle=:dot)
        vline!(p, [tinf_k], color=:gray, lw=0.8, linestyle=:dash, label="")
        push!(panels_4b, p)
    end
end

for N0_n in N0_waarden
    lbl = "N0=$(Int(round(N0_n))), MOI=2.0, t_inf=2h"
    p = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
        legend=:outertopright,
        tickfontsize=9, guidefontsize=10, legendfontsize=8,
        bottom_margin=6Plots.mm, left_margin=8Plots.mm,
        title=lbl, titlefontsize=9)

    s3 = Model3.run(make_p3(N0_n, 2.0, 2.0))
    s4 = Model4.run(make_p4(N0_n, 2.0, 2.0))
    s5 = Model5.run(make_p5(N0_n, 2.0, 2.0))

    SciMLBase.successful_retcode(s3) && plot!(p, s3.t,
        [s3.u[i][Model3.Benzind] for i in eachindex(s3.u)],
        label="M3", color=kleur_m3, lw=2)
    SciMLBase.successful_retcode(s4) && plot!(p, s4.t,
        [s4.u[i][Model4.Benzind] for i in eachindex(s4.u)],
        label="M4", color=kleur_m4, lw=2, linestyle=:dash)
    SciMLBase.successful_retcode(s5) && plot!(p, s5.t,
        [s5.u[i][Model5.Benzind] for i in eachindex(s5.u)],
        label="M5", color=kleur_m5, lw=2, linestyle=:dot)
    vline!(p, [2.0], color=:gray, lw=0.8, linestyle=:dash, label="")
    push!(panels_4b, p)
end

fig4b = plot(panels_4b..., layout=(4,3), size=(1500,1600), margin=5Plots.mm)
savefig(fig4b, "h4b_benzonase_drie_modellen.png")
println("  Saved: h4b_benzonase_drie_modellen.png")

# ============================================================
#  4e. Modelvergelijking als functie van parameters
# ============================================================
println("\n=== 4e: Modelvergelijking als functie van parameters ===")

moi_sweep_e    = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
tinf_sweep_e   = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]
N0_sweep_e     = [1e6, 1e7, 1e8, 5e8, 1e9, 2e9, 5e9]

benz3_moi = Float64[]; benz4_moi = Float64[]; benz5_moi = Float64[]
for moi in moi_sweep_e
    s3 = Model3.run(make_p3(1e9, 2.0, moi))
    s4 = Model4.run(make_p4(1e9, 2.0, moi))
    s5 = Model5.run(make_p5(1e9, 2.0, moi))
    push!(benz3_moi, SciMLBase.successful_retcode(s3) ? maximum([s3.u[i][Model3.Benzind] for i in eachindex(s3.u)]) : NaN)
    push!(benz4_moi, SciMLBase.successful_retcode(s4) ? maximum([s4.u[i][Model4.Benzind] for i in eachindex(s4.u)]) : NaN)
    push!(benz5_moi, SciMLBase.successful_retcode(s5) ? maximum([s5.u[i][Model5.Benzind] for i in eachindex(s5.u)]) : NaN)
end

benz3_tinf = Float64[]; benz4_tinf = Float64[]; benz5_tinf = Float64[]
for tinf in tinf_sweep_e
    s3 = Model3.run(make_p3(1e9, tinf, 2.0))
    s4 = Model4.run(make_p4(1e9, tinf, 2.0))
    s5 = Model5.run(make_p5(1e9, tinf, 2.0))
    push!(benz3_tinf, SciMLBase.successful_retcode(s3) ? maximum([s3.u[i][Model3.Benzind] for i in eachindex(s3.u)]) : NaN)
    push!(benz4_tinf, SciMLBase.successful_retcode(s4) ? maximum([s4.u[i][Model4.Benzind] for i in eachindex(s4.u)]) : NaN)
    push!(benz5_tinf, SciMLBase.successful_retcode(s5) ? maximum([s5.u[i][Model5.Benzind] for i in eachindex(s5.u)]) : NaN)
end

benz3_N0 = Float64[]; benz4_N0 = Float64[]; benz5_N0 = Float64[]
for N0_ in N0_sweep_e
    s3 = Model3.run(make_p3(N0_, 2.0, 2.0))
    s4 = Model4.run(make_p4(N0_, 2.0, 2.0))
    s5 = Model5.run(make_p5(N0_, 2.0, 2.0))
    push!(benz3_N0, SciMLBase.successful_retcode(s3) ? maximum([s3.u[i][Model3.Benzind] for i in eachindex(s3.u)]) : NaN)
    push!(benz4_N0, SciMLBase.successful_retcode(s4) ? maximum([s4.u[i][Model4.Benzind] for i in eachindex(s4.u)]) : NaN)
    push!(benz5_N0, SciMLBase.successful_retcode(s5) ? maximum([s5.u[i][Model5.Benzind] for i in eachindex(s5.u)]) : NaN)
end

println("\n  Variatie tussen modellen (max relatief verschil):")
@printf("  MOI sweep:   M3 range=[%.4e, %.4e] | M4=[%.4e, %.4e] | M5=[%.4e, %.4e]\n",
    minimum(filter(!isnan, benz3_moi)), maximum(filter(!isnan, benz3_moi)),
    minimum(filter(!isnan, benz4_moi)), maximum(filter(!isnan, benz4_moi)),
    minimum(filter(!isnan, benz5_moi)), maximum(filter(!isnan, benz5_moi)))
@printf("  t_inf sweep: M3 range=[%.4e, %.4e] | M4=[%.4e, %.4e] | M5=[%.4e, %.4e]\n",
    minimum(filter(!isnan, benz3_tinf)), maximum(filter(!isnan, benz3_tinf)),
    minimum(filter(!isnan, benz4_tinf)), maximum(filter(!isnan, benz4_tinf)),
    minimum(filter(!isnan, benz5_tinf)), maximum(filter(!isnan, benz5_tinf)))
@printf("  N0 sweep:    M3 range=[%.4e, %.4e] | M4=[%.4e, %.4e] | M5=[%.4e, %.4e]\n",
    minimum(filter(!isnan, benz3_N0)), maximum(filter(!isnan, benz3_N0)),
    minimum(filter(!isnan, benz4_N0)), maximum(filter(!isnan, benz4_N0)),
    minimum(filter(!isnan, benz5_N0)), maximum(filter(!isnan, benz5_N0)))

p4e_moi = plot(xlabel="Initial MOI [-]", ylabel="Max Benzonase [mmol L⁻¹]",
    xscale=:log10, legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
plot!(p4e_moi, moi_sweep_e, benz3_moi, marker=:circle,  lw=2, color=:steelblue,  label="Model 3")
plot!(p4e_moi, moi_sweep_e, benz4_moi, marker=:square,  lw=2, color=:darkorange, label="Model 4", linestyle=:dash)
plot!(p4e_moi, moi_sweep_e, benz5_moi, marker=:diamond, lw=2, color=:darkgreen,  label="Model 5", linestyle=:dot)

p4e_tinf = plot(xlabel="t_inf [h]", ylabel="Max Benzonase [mmol L⁻¹]",
    legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
plot!(p4e_tinf, tinf_sweep_e, benz3_tinf, marker=:circle,  lw=2, color=:steelblue,  label="Model 3")
plot!(p4e_tinf, tinf_sweep_e, benz4_tinf, marker=:square,  lw=2, color=:darkorange, label="Model 4", linestyle=:dash)
plot!(p4e_tinf, tinf_sweep_e, benz5_tinf, marker=:diamond, lw=2, color=:darkgreen,  label="Model 5", linestyle=:dot)

p4e_N0 = plot(xlabel="N₀ [cells L⁻¹]", ylabel="Max Benzonase [mmol L⁻¹]",
    xscale=:log10, legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
plot!(p4e_N0, N0_sweep_e, benz3_N0, marker=:circle,  lw=2, color=:steelblue,  label="Model 3")
plot!(p4e_N0, N0_sweep_e, benz4_N0, marker=:square,  lw=2, color=:darkorange, label="Model 4", linestyle=:dash)
plot!(p4e_N0, N0_sweep_e, benz5_N0, marker=:diamond, lw=2, color=:darkgreen,  label="Model 5", linestyle=:dot)

fig4e = plot(p4e_moi, p4e_tinf, p4e_N0,
    layout=(1,3), size=(1400,430), margin=6Plots.mm)
savefig(fig4e, "h4e_modelvergelijking_parameters.png")
println("  Saved: h4e_modelvergelijking_parameters.png")

# ============================================================
#  4c. Phage population dynamics
# ============================================================
println("\n=== 4c: Phage population dynamics ===")

fig4c = plot(xlabel="t [h]", ylabel="Free phages L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:outertopright,
    tickfontsize=11, guidefontsize=13, legendfontsize=10,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(fig4c, t3, Pf3, label="Model 3 (fixed b)",  color=:steelblue,  lw=2)
plot!(fig4c, t4, Pf4, label="Model 4 (linear b)", color=:darkorange, lw=2, linestyle=:dash)
plot!(fig4c, t5, Pf5, label="Model 5 (FBA b)",     color=:darkgreen,  lw=2, linestyle=:dot)
savefig(fig4c, "h4c_fagen_drie_modellen.png")
println("  Saved: h4c_fagen_drie_modellen.png")

# ============================================================
#  4d. Summary table
# ============================================================
println("\n=== 4d: Summary table ===")

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

println("\n  Model   | Max Benz [mmol/L] | P/N end   | Score (0.8B-0.2PN)")
for (naam, benz, pn) in [("Model 3", maximum(B3), pn3), ("Model 4", maximum(B4), pn4), ("Model 5", maximum(B5), pn5)]
    score = 0.8*(benz/benz_ref) - 0.2*(isinf(pn) ? 1.0 : pn/pn_ref)
    @printf("  %-7s | %-17.5f | %-9.2e | %.4f\n", naam, benz, pn, score)
end

fig4d = bar(["Model 3" "Model 4" "Model 5"],
    [maximum(B3) maximum(B4) maximum(B5)],
    ylabel="Max Benzonase [mmol L⁻¹]",
    color=[:steelblue :darkorange :darkgreen],
    tickfontsize=11, guidefontsize=13,
    legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
savefig(fig4d, "h4d_samenvatting_drie_modellen.png")
println("  Saved: h4d_samenvatting_drie_modellen.png")

println("\n=== Chapter 4 complete ===")