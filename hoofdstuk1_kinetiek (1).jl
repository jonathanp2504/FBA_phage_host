# ============================================================
#  HOOFDSTUK 1: Kinetische validatie Model 1
#  Drie analyses:
#   1a. Één-stap groeicurve: latente periode en burst size
#   1b. MOI-sweep: lysis-lysogenie beslissing
#   1c. Gevoeligheid voor alfa_ads en tau
# ============================================================

include("./Bin_model_1/parameters.jl")
include("./Bin_model_1/FBA.jl")
include("./Bin_model_1/dFBA.jl")

using Plots, Statistics, SciMLBase
using COBREXA, AbstractFBCModels
import SBMLFBCModels

model_path = joinpath(@__DIR__, "iJO1366.xml")
model      = loadFBAmodel(model_path)

alpha_syn      = 2.0;  beta_deg = 0.5
K_s            = [0.0278, 0.0146, 0.0543, 0.0833]
tau            = 1.0;  b = 170.0;  alfa_ads = 1e-10
E_coli_cellDW  = 2.8e-13
p_pref         = [0.8925, 0.08925, 0.008925, 0.008925]
V_max          = [12.7, 3.75, 0.0, 4.0]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
essentials_ids = ["R_EX_o2_e","R_EX_nh4_e","R_EX_pi_e","R_EX_so4_e",
                  "R_EX_k_e","R_EX_mg2_e","R_EX_ca2_e","R_EX_cl_e",
                  "R_EX_fe2_e","R_EX_fe3_e","R_EX_mn2_e","R_EX_zn2_e",
                  "R_EX_cu2_e","R_EX_cobalt2_e","R_EX_mobd_e","R_EX_thi_e",
                  "R_EX_ni2_e","R_EX_sel_e","R_EX_slnt_e","R_EX_tungs_e"]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release      = 5.66e-12;  duration = 15.0
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)
fbaCache       = buildFbaCache(model, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")

function make_p1(N0, t_inf, moi; tau_=tau, alfa_=alfa_ads)
    Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau_, b, alfa_, E_coli_cellDW, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        fbaCache, 0.0, zeros(4), t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0, mu_max_vector, e_max_vector, 0.0, 0.0)
end

# ============================================================
#  1a. Één-stap groeicurve validatie (lage MOI = lytisch regime)
# ============================================================
println("=== 1a: Groeicurve validatie ===")
p_ref = make_p1(1e9, 2.0, 0.001)
sol   = run(p_ref)

t   = sol.t
S   = [sol.u[i][Sind_S] for i in eachindex(sol.u)]
I   = [sol.u[i][Sind_I] for i in eachindex(sol.u)]
L   = [sol.u[i][Sind_L] for i in eachindex(sol.u)]
P   = [sol.u[i][Pfind]  for i in eachindex(sol.u)]
MOI = [min(sol.u[i][Pfind] / max(sol.u[i][Sind_S], 1.0), 10.0)
       for i in eachindex(sol.u)]
X   = S .+ I .+ L

i_max_I = argmax(I)
t_lysis = t[i_max_I]
println("  Lysis-piek bij t = $(round(t_lysis,digits=2)) h (verwacht ≈ $(2.0+tau) h)")
println("  Ingestelde burst size b = $b | tau = $tau h | alfa_ads = $alfa_ads L/cel/h")

pa = plot(t, [S I L X],
    label=["Vatbaar S" "Lytisch I" "Lysogeen L" "Totaal X"],
    color=[:blue :red :green :black], lw=2,
    yscale=:log10, ylims=(1,:auto),
    title="1a: Populatiedynamica (MOI=0.001, lytisch regime)",
    xlabel="t [h]", ylabel="cellen/L")
pb = plot(t, max.(P, 1.0), label="Vrije fagen P", color=:darkred, lw=2,
    yscale=:log10, ylims=(1,:auto),
    title="1a: Faagpopulatie", xlabel="t [h]", ylabel="fagen/L")
pc_moi = plot(t, MOI .+ 1e-6, label="Effectieve MOI", color=:purple, lw=2,
    yscale=:log10,
    title="1a: MOI tijdsverloop", xlabel="t [h]", ylabel="MOI [-]")
fig1a = plot(pa, pb, pc_moi, layout=(1,3), size=(1200,380))
savefig(fig1a, "h1a_groeicurve_validatie.png")
println("  Figuur opgeslagen: h1a_groeicurve_validatie.png")

# ============================================================
#  1b. MOI-sweep: lysis-lysogenie beslissing
# ============================================================
println("\n=== 1b: MOI-sweep ===")
moi_sweep  = [0.0001, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
frac_lys   = Float64[]
frac_lyt   = Float64[]
final_P_sw = Float64[]

for moi in moi_sweep
    sol_sw = run(make_p1(1e9, 2.0, moi))
    if SciMLBase.successful_retcode(sol_sw)
        # Fracties meten op tijdstip van maximale infectie (niet einde)
        infected = [sol_sw.u[i][Sind_I] + sol_sw.u[i][Sind_L]
                    for i in eachindex(sol_sw.u)]
        i_max = argmax(infected)
        I_f   = sol_sw.u[i_max][Sind_I]
        L_f   = sol_sw.u[i_max][Sind_L]
        inf_f = I_f + L_f
        push!(frac_lys,   inf_f > 1.0 ? L_f/inf_f : 0.0)
        push!(frac_lyt,   inf_f > 1.0 ? I_f/inf_f : 0.0)
        push!(final_P_sw, max(sol_sw.u[end][Pfind], 1.0))
        println("  MOI=$moi: frac_lys=$(round(frac_lys[end],digits=3)) | ",
                "frac_lyt=$(round(frac_lyt[end],digits=3)) | ",
                "P_eind=$(round(final_P_sw[end], sigdigits=2))")
    else
        println("  MOI=$moi: simulatie gefaald")
        push!(frac_lys,   0.0)
        push!(frac_lyt,   0.0)
        push!(final_P_sw, 1.0)
    end
end

pd = plot(moi_sweep, [frac_lys frac_lyt],
    label=["Fractie lysogeen" "Fractie lytisch"],
    color=[:green :red], marker=:circle, lw=2,
    xscale=:log10, ylims=(0,1.05),
    xlabel="MOI", ylabel="Fractie geïnfecteerde cellen",
    title="1b: Lysis-Lysogenie beslissing als functie van MOI")
pe = plot(moi_sweep, final_P_sw,
    label="Vrije fagen (eindwaarde)", color=:darkred,
    marker=:circle, lw=2, xscale=:log10, yscale=:log10,
    xlabel="MOI", ylabel="fagen/L",
    title="1b: Eindconcentratie fagen vs MOI")
fig1b = plot(pd, pe, layout=(1,2), size=(900,380))
savefig(fig1b, "h1b_moi_sweep.png")
println("  Figuur opgeslagen: h1b_moi_sweep.png")

# ============================================================
#  1c. Gevoeligheidsanalyse: alfa_ads en tau
# ============================================================
println("\n=== 1c: Gevoeligheidsanalyse ===")

alfa_waarden = [1e-11, 1e-10, 5e-10, 1e-9]
tau_waarden  = [0.5, 1.0, 1.5, 2.0]

max_P_alfa = Float64[]
for alfa in alfa_waarden
    sol_a = run(make_p1(1e9, 2.0, 0.001; alfa_=alfa))
    push!(max_P_alfa, SciMLBase.successful_retcode(sol_a) ?
        maximum([sol_a.u[i][Pfind] for i in eachindex(sol_a.u)]) : 0.0)
    println("  alfa_ads=$alfa: max_P=$(round(max_P_alfa[end], sigdigits=2))")
end

t_piek_tau = Float64[]
for tau_ in tau_waarden
    sol_t = run(make_p1(1e9, 2.0, 0.001; tau_=tau_))
    if SciMLBase.successful_retcode(sol_t)
        I_t = [sol_t.u[i][Sind_I] for i in eachindex(sol_t.u)]
        push!(t_piek_tau, sol_t.t[argmax(I_t)])
    else
        push!(t_piek_tau, NaN)
    end
    println("  tau=$(tau_) h: lysis-piek bij t=$(round(t_piek_tau[end],digits=2)) h")
end

pf = bar(string.(alfa_waarden), max_P_alfa,
    title="1c: Max faagconcentratie vs adsorptieconstante",
    xlabel="alfa_ads [L/cel/h]", ylabel="Max fagen/L",
    color=:steelblue, legend=false)
pg = bar(string.(tau_waarden), t_piek_tau,
    title="1c: Lysis-piek tijdstip vs latente periode",
    xlabel="tau [h]", ylabel="t_piek [h]",
    color=:darkorange, legend=false)
fig1c = plot(pf, pg, layout=(1,2), size=(900,380))
savefig(fig1c, "h1c_gevoeligheidsanalyse.png")
println("  Figuur opgeslagen: h1c_gevoeligheidsanalyse.png")

println("\n=== Hoofdstuk 1 voltooid ===")
