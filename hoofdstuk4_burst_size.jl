# ============================================================
#  HOOFDSTUK 4: Burst size — vast, lineair of metabolisch?
#  Model 3 (vast b) vs Model 4 (lineair b(µ)) vs Model 5 (FBA b)
#
#  Analyses:
#   4a. Burst size tijdsverloop TOT takeover + substraatverloop
#   4b. Benzonase vergelijking: 3×3 condities (MOI × t_inf)
#       om te bevestigen dat Model 3/4/5 nauwelijks verschilt qua productie
#   4c. Faagpopulatiedynamica
#   4d. Samenvattingstabel
# ============================================================
include("Model3.jl")
include("Model4.jl")
include("Model5.jl")

using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels

model_path     = joinpath(@__DIR__, "iJO1366.xml")
alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0061, 9.4e-4, 0.0543, 8.33]
tau            = 1.32
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
N0=1e9; t_inf=2.0; moi=2.0; f_prod=0.0015

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

function make_p3(N0_, t_inf_, moi_, f_=f_prod)
    Model3.Parameters(duration, N0_, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, 170.0, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba3, lysogenFba3,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf_, moi_*N0_,
        "R_BENZ_prod", 0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_)
end

function make_p4(N0_, t_inf_, moi_, f_=f_prod)
    Model4.Parameters(duration, N0_, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, 50.0, 228.6, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba4, lysogenFba4,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf_, moi_*N0_,
        "R_BENZ_prod", 0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_)
end

function make_p5(N0_, t_inf_, moi_, f_=f_prod)
    Model5.Parameters(duration, N0_, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba5, lysogenFba5, lyticFba5,
        0.0, zeros(4), 0.0, zeros(4), 0.0, 0.0,
        7.92e-8, 6.48, 3.02, 0.01, t_inf_, moi_*N0_,
        "R_BENZ_prod", "R_PHAGE_prod",
        0.05, 0.1, 0.0, mu_max_vector, e_max_vector, f_)
end

# ============================================================
#  Referentiesimulaties bij standaardcondities
# ============================================================
println("=== Referentiesimulaties (N0=1e9, t_inf=2h, MOI=2) ===")

p3 = make_p3(N0, t_inf, moi)
p4 = make_p4(N0, t_inf, moi)
p5 = make_p5(N0, t_inf, moi)

sol3 = Model3.run(p3)
sol4 = Model4.run(p4)
sol5 = Model5.run(p5)

t3 = sol3.t; B3 = [sol3.u[i][Model3.Benzind] for i in eachindex(sol3.u)]
              Pf3= [max(sol3.u[i][Model3.Pfind], 1.0) for i in eachindex(sol3.u)]
              Lys3=[sol3.u[i][Model3.lind] for i in eachindex(sol3.u)]
              N3  =[sol3.u[i][Model3.Nind] for i in eachindex(sol3.u)]
              S3  =[sol3.u[i][Model3.Sind[2]] for i in eachindex(sol3.u)]  # maltose

t4 = sol4.t; B4 = [sol4.u[i][Model4.Benzind] for i in eachindex(sol4.u)]
              Pf4= [max(sol4.u[i][Model4.Pfind], 1.0) for i in eachindex(sol4.u)]
              N4  =[sol4.u[i][Model4.Nind] for i in eachindex(sol4.u)]

t5 = sol5.t; B5 = [sol5.u[i][Model5.Benzind] for i in eachindex(sol5.u)]
              Pf5= [max(sol5.u[i][Model5.Pfind], 1.0) for i in eachindex(sol5.u)]

# Burst size tijdsreeksen
b_tijds_m3 = fill(170.0, length(t3))
# Model 4: reconstrueer via numerieke differentiatie op N4
mu_N_recon = zeros(length(t4))
for i in 2:length(t4)-1
    dt = t4[i+1] - t4[i-1]
    if N4[i] > 1.0 && N4[i-1] > 1.0 && dt > 0
        mu_N_recon[i] = max(0.0, log(max(N4[i+1],1.0)/max(N4[i-1],1.0)) / dt)
    end
end
mu_N_recon[1] = mu_N_recon[2]; mu_N_recon[end] = mu_N_recon[end-1]
b_tijds_m4 = 50.0 .+ 228.6 .* mu_N_recon
b_tijds_m5 = [sol5.u[i][Model5.Bind] for i in eachindex(sol5.u)]

println("  Max Benzonase M3: $(round(maximum(B3),sigdigits=3)) | M4: $(round(maximum(B4),sigdigits=3)) | M5: $(round(maximum(B5),sigdigits=3))")

# ============================================================
#  4a. Burst size tijdsverloop TOT takeover + substraatverloop
#
#  "Takeover" = tijdstip waarop naïeve cellen < 5% van totale biomassa
#  Na takeover heeft burst size geen biologische relevantie meer omdat
#  bijna alle cellen al lysogeen zijn en niet meer geïnfecteerd worden.
# ============================================================
println("\n=== 4a: Burst size tijdsverloop (tot takeover) ===")

# Bepaal takeover-tijdstip: N < 5% van totale biomassa
function find_takeover(sol, Nind_, Lind_, lind_)
    for i in eachindex(sol.u)
        N_  = sol.u[i][Nind_]
        L_  = sol.u[i][Lind_]
        l_  = sol.u[i][lind_]
        X_  = N_ + L_ + l_
        if X_ > 1e4 && N_ / X_ < 0.05
            return sol.t[i]
        end
    end
    return sol.t[end]
end

t_takeover3 = find_takeover(sol3, Model3.Nind, Model3.Lind, Model3.lind)
t_takeover4 = find_takeover(sol4, Model4.Nind, Model4.Lind, Model4.lind)
t_takeover5 = find_takeover(sol5, Model5.Nind, Model5.Lind, Model5.lind)
t_takeover  = minimum([t_takeover3, t_takeover4, t_takeover5])
println("  Takeover tijdstippen: M3=$(round(t_takeover3,digits=2)) | M4=$(round(t_takeover4,digits=2)) | M5=$(round(t_takeover5,digits=2)) h")

# Filter tot takeover
idx3 = findall(t3 .<= t_takeover)
idx4 = findall(t4 .<= t_takeover)
idx5 = findall(t5 .<= t_takeover)

# Panel A: burst size tot takeover
p4a_burst = plot(xlabel="t [h]", ylabel="Burst size b [fagen cel⁻¹]",
    ylims=(0, max(maximum(b_tijds_m4[idx4])*1.15, 250)),
    legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p4a_burst, t3[idx3], b_tijds_m3[idx3],
    label="Model 3 (vast b=170)", color=:steelblue, lw=2, linestyle=:dash)
plot!(p4a_burst, t4[idx4], b_tijds_m4[idx4],
    label="Model 4 (lineair b(µ))", color=:darkorange, lw=2)
plot!(p4a_burst, t5[idx5], b_tijds_m5[idx5],
    label="Model 5 (FBA b)", color=:darkgreen, lw=2, linestyle=:dot)
vline!(p4a_burst, [t_inf], color=:gray, lw=1, linestyle=:dash, label="t_inf")
vline!(p4a_burst, [t_takeover], color=:black, lw=1.5, linestyle=:dot, label="Takeover")

# Panel B: substraatverloop (maltose als proxy voor LamB-beschikbaarheid)
S_malt3 = [sol3.u[i][Model3.Sind[2]] for i in eachindex(sol3.u)]  # maltose
S_glc3  = [sol3.u[i][Model3.Sind[1]] for i in eachindex(sol3.u)]  # glucose
S_glyc3 = [sol3.u[i][Model3.Sind[3]] for i in eachindex(sol3.u)]  # glycerol
S_ac3   = [sol3.u[i][Model3.Sind[4]] for i in eachindex(sol3.u)]  # acetaat

p4a_subs = plot(t3[idx3], S_glc3[idx3],
    label="Glucose", color=:steelblue, lw=2,
    xlabel="t [h]", ylabel="Substraatconcentratie [mmol L⁻¹]",
    legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(p4a_subs, t3[idx3], S_malt3[idx3],
    label="Maltose", color=:darkorange, lw=2, linestyle=:dash)
plot!(p4a_subs, t3[idx3], S_glyc3[idx3],
    label="Glycerol", color=:green, lw=2, linestyle=:dot)
plot!(p4a_subs, t3[idx3], S_ac3[idx3],
    label="Acetaat", color=:purple, lw=1.5, linestyle=:dashdot)
vline!(p4a_subs, [t_inf], color=:gray, lw=1, linestyle=:dash, label="t_inf")

fig4a = plot(p4a_burst, p4a_subs, layout=(1,2), size=(1100,430), margin=6Plots.mm)
savefig(fig4a, "h4a_burst_size_tijdsverloop.png")
println("  Figuur opgeslagen: h4a_burst_size_tijdsverloop.png")

# ============================================================
#  4b. Benzonase vergelijking: 3×3 condities (MOI × t_inf)
#  beginbiomassa = 1e9 in alle gevallen
#  MOI  = [0.1, 0.5, 2.0]
#  t_inf = [1.0, 2.0, 5.0] h
# ============================================================
println("\n=== 4b: Benzonase vergelijking 3×3 condities ===")

moi_waarden  = [0.1, 0.5, 2.0]
tinf_waarden = [1.0, 2.0, 5.0]
kleuren_4b   = [:steelblue :darkorange :darkgreen]
stijlen_4b   = [:solid :dash :dot]

fig4b_panels = []

for (j, moi_j) in enumerate(moi_waarden)
    p_benz = plot(xlabel="t [h]", ylabel="Benzonase [mmol L⁻¹]",
        legend= j==1 ? :topleft : false,
        bottom_margin=6Plots.mm, left_margin=10Plots.mm,
        title="MOI=$(moi_j)")

    for (k, tinf_k) in enumerate(tinf_waarden)
        sol3_jk = Model3.run(make_p3(1e9, tinf_k, moi_j))
        sol4_jk = Model4.run(make_p4(1e9, tinf_k, moi_j))
        sol5_jk = Model5.run(make_p5(1e9, tinf_k, moi_j))

        lbl3 = j==1 && k==1 ? "M3 t_inf=$(tinf_k)h" : (j==1 ? "M3 t_inf=$(tinf_k)h" : "")
        lbl4 = j==1 && k==1 ? "M4 t_inf=$(tinf_k)h" : (j==1 ? "M4 t_inf=$(tinf_k)h" : "")
        lbl5 = j==1 && k==1 ? "M5 t_inf=$(tinf_k)h" : (j==1 ? "M5 t_inf=$(tinf_k)h" : "")

        if SciMLBase.successful_retcode(sol3_jk)
            B3_jk = [sol3_jk.u[i][Model3.Benzind] for i in eachindex(sol3_jk.u)]
            plot!(p_benz, sol3_jk.t, B3_jk,
                color=kleuren_4b[k], lw=2,   linestyle=:solid, label=lbl3)
        end
        if SciMLBase.successful_retcode(sol4_jk)
            B4_jk = [sol4_jk.u[i][Model4.Benzind] for i in eachindex(sol4_jk.u)]
            plot!(p_benz, sol4_jk.t, B4_jk,
                color=kleuren_4b[k], lw=1.5, linestyle=:dash,  label=lbl4)
        end
        if SciMLBase.successful_retcode(sol5_jk)
            B5_jk = [sol5_jk.u[i][Model5.Benzind] for i in eachindex(sol5_jk.u)]
            plot!(p_benz, sol5_jk.t, B5_jk,
                color=kleuren_4b[k], lw=1.5, linestyle=:dot,   label=lbl5)
        end

        # Voeg t_inf annotatie toe
        vline!(p_benz, [tinf_k], color=kleuren_4b[k], lw=0.8, linestyle=:dashdot, label="")
    end
    push!(fig4b_panels, p_benz)
end

# Legenda uitleg: kleur = t_inf, lijnstijl = model
# solid=M3, dash=M4, dot=M5; blauw=1h, oranje=2h, groen=5h
fig4b = plot(fig4b_panels..., layout=(1,3), size=(1400,440), margin=6Plots.mm)
savefig(fig4b, "h4b_benzonase_drie_modellen.png")
println("  Figuur opgeslagen: h4b_benzonase_drie_modellen.png")
println("  Legenda: kleur = t_inf (blauw=1h, oranje=2h, groen=5h)")
println("           lijnstijl = model (solid=M3, dash=M4, dot=M5)")

# ============================================================
#  4c. Faagpopulatiedynamica (referentiecondities)
# ============================================================
println("\n=== 4c: Faagpopulatiedynamica ===")

fig4c = plot(xlabel="t [h]", ylabel="Vrije fagen L⁻¹",
    yscale=:log10, ylims=(1,:auto), legend=:topright,
    bottom_margin=6Plots.mm, left_margin=10Plots.mm)
plot!(fig4c, t3, Pf3, label="Model 3 (vast b)",    color=:steelblue,  lw=2)
plot!(fig4c, t4, Pf4, label="Model 4 (lineair b)", color=:darkorange, lw=2, linestyle=:dash)
plot!(fig4c, t5, Pf5, label="Model 5 (FBA b)",     color=:darkgreen,  lw=2, linestyle=:dot)
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

results = [("Model 3", maximum(B3), pn3),
           ("Model 4", maximum(B4), pn4),
           ("Model 5", maximum(B5), pn5)]

println("\n  Model   | Max Benz [mmol/L] | P/N eind  | Score (0.8B-0.2PN)")
println("  --------|-------------------|-----------|--------------------")
for (naam, benz, pn) in results
    score = 0.8*(benz/benz_ref) - 0.2*(isinf(pn) ? 1.0 : pn/pn_ref)
    @printf("  %-7s | %-17.5f | %-9.2e | %.4f\n", naam, benz, pn, score)
end

fig4d = bar(["Model 3" "Model 4" "Model 5"],
    [maximum(B3) maximum(B4) maximum(B5)],
    ylabel="Max Benzonase [mmol L⁻¹]",
    color=[:steelblue :darkorange :darkgreen],
    legend=false,
    bottom_margin=8Plots.mm, left_margin=10Plots.mm)
savefig(fig4d, "h4d_samenvatting_drie_modellen.png")
println("  Figuur opgeslagen: h4d_samenvatting_drie_modellen.png")

println("\n=== Hoofdstuk 4 voltooid ===")