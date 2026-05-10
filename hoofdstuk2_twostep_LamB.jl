# ============================================================
#  HOOFDSTUK 2: Two-step adsorptie en de rol van LamB
#  Vergelijking Model 1 vs Model 2 + glucose-maltose puls
# ============================================================
include("Model1.jl")
include("Model2.jl")
using Plots, Statistics, Printf, SciMLBase
using COBREXA, AbstractFBCModels, DelayDiffEq, OrdinaryDiffEq
import SBMLFBCModels, OrdinaryDiffEqCore



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
h_release      = 5.66e-12;  duration = 15.0
mu_max_vector  = [1.33, 1.26, 1.10, 0.29]
e_max_vector   = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)

# FBA caches
baseModel2  = Model2.loadFBAmodel(model_path)
all_ex_ids  = [id for id in keys(baseModel2.reactions) if startswith(id, "R_EX_")]
naiveFba2   = Model2.buildFbaCache(baseModel2, exchange_ids, "R_BIOMASS_Ec_iJO1366_core_53p95M")
lysogenFba2 = Model2.buildFbaCache(Model2.loadFBAmodel(model_path), exchange_ids,
                                    "R_BIOMASS_Ec_iJO1366_core_53p95M")
fbaCache1   = Model1.buildFbaCache(Model1.loadFBAmodel(model_path), exchange_ids,
                                    "R_BIOMASS_Ec_iJO1366_core_53p95M")

function make_p2(N0, t_inf, moi)
    Model2.Parameters(duration, N0, alpha_syn, beta_deg, K_s, V_max, p_pref,
        tau, b, 1e-12, MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        naiveFba2, lysogenFba2,
        0.0, zeros(4), 0.0, zeros(4), 0.0,
        1e-10, 10.0, 5.0, 0.01,
        t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0,
        mu_max_vector, e_max_vector, 0.0, 0.0)
end

function make_p1(N0, t_inf, moi)
    Model1.Parameters(duration, N0, alpha_syn, beta_deg, K_s,
        [12.7,3.75,0.0,4.0], p_pref, tau, b, 1e-10, 2.8e-13,
        MW_values, h_release,
        "R_BIOMASS_Ec_iJO1366_core_53p95M", exchange_ids, all_ex_ids, essentials_ids,
        fbaCache1, 0.0, zeros(4), t_inf, moi*N0,
        "R_BENZ_prod", 0.0, 0.0, 0.0, mu_max_vector, e_max_vector, 0.0, 0.0)
end

# ============================================================
#  2a. Populatiedynamica Model 2 bij 4 MOI-waarden
# ============================================================
println("=== 2a: Model 2 populatiedynamica bij verschillende MOI ===")

moi_compare = [0.001, 0.1, 1.0, 5.0]
panels = []

for moi in moi_compare
    sol2 = Model2.run(make_p2(1e9, 2.0, moi))
    t2   = sol2.t
    N2   = [sol2.u[i][Model2.Nind]  for i in eachindex(sol2.u)]
    Lyt2 = [sol2.u[i][Model2.Lind]  for i in eachindex(sol2.u)]
    Lys2 = [sol2.u[i][Model2.lind]  for i in eachindex(sol2.u)]
    Pf2  = [max(sol2.u[i][Model2.Pfind], 1.0) for i in eachindex(sol2.u)]
    Pa2  = [sol2.u[i][Model2.Paind] for i in eachindex(sol2.u)]

    p_panel = plot(t2, [N2 Lyt2 Lys2 Pf2],
        label=["Naïef N" "Lytisch L" "Lysogeen l" "Vrije fagen Pf"],
        color=[:blue :red :green :darkred],
        lw=2, yscale=:log10, ylims=(1,:auto),
        title="MOI = $moi", xlabel="t [h]", ylabel="cellen of fagen/L",
        legend = moi==moi_compare[1] ? :topright : false)
    plot!(p_panel, t2, max.(Pa2, 1.0), label="Gehecht Pa",
        color=:purple, lw=1.5, linestyle=:dot,
        legend = moi==moi_compare[1] ? :topright : false)
    push!(panels, p_panel)
end

fig2a = plot(panels..., layout=(2,2), size=(1000,700),
    plot_title="2a: Model 2 populatiedynamica bij verschillende MOI")
savefig(fig2a, "h2a_model2_moi_vergelijking.png")
println("  Figuur opgeslagen: h2a_model2_moi_vergelijking.png")

# ============================================================
#  2b. Lysogene fractie: Model 1 vs Model 2
# ============================================================
println("\n=== 2b: Lysogene fractie Model 1 vs Model 2 ===")

moi_sweep = [0.0001, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
fl_m1 = Float64[]; fl_m2 = Float64[]; valid_moi = Float64[]

for moi in moi_sweep
    sol1 = Model1.run(make_p1(1e9, 2.0, moi))
    sol2 = Model2.run(make_p2(1e9, 2.0, moi))

    if SciMLBase.successful_retcode(sol1) && SciMLBase.successful_retcode(sol2)
        # Model 1: op tijdstip van maximale infectie
        infected1 = [sol1.u[i][Model1.Sind_I] + sol1.u[i][Model1.Sind_L]
                     for i in eachindex(sol1.u)]
        i_max1 = argmax(infected1)
        I1 = sol1.u[i_max1][Model1.Sind_I]; L1 = sol1.u[i_max1][Model1.Sind_L]
        inf1 = I1 + L1
        push!(fl_m1, inf1 > 1.0 ? L1/inf1 : 0.0)

        # Model 2: op tijdstip van maximale infectie
        infected2 = [sol2.u[i][Model2.Lind] + sol2.u[i][Model2.lind]
                     for i in eachindex(sol2.u)]
        i_max2 = argmax(infected2)
        Lyt2 = sol2.u[i_max2][Model2.Lind]; Lys2 = sol2.u[i_max2][Model2.lind]
        inf2 = Lyt2 + Lys2
        push!(fl_m2, inf2 > 1.0 ? Lys2/inf2 : 0.0)

        push!(valid_moi, moi)
        println("  MOI=$moi: M1=$(round(fl_m1[end],digits=3)) | M2=$(round(fl_m2[end],digits=3))")
    end
end

fig2b = plot(valid_moi, [fl_m1 fl_m2],
    label=["Model 1 (one-step)" "Model 2 (two-step LamB)"],
    marker=[:circle :square], lw=2, color=[:steelblue :darkorange],
    xscale=:log10, ylims=(0,1.05),
    xlabel="MOI", ylabel="Fractie lysogeen",
    title="2b: Lysis-Lysogenie: Model 1 vs Model 2",
    legend=:bottomright, size=(650,420))
savefig(fig2b, "h2b_lysogenie_model1_vs_model2.png")
println("  Figuur opgeslagen: h2b_lysogenie_model1_vs_model2.png")

println("\n  MOI   | Model1 | Model2 | Verschil")
for idx in eachindex(valid_moi)
    @printf("  %-5.4f | %-6.3f | %-6.3f | %.3f\n",
        valid_moi[idx], fl_m1[idx], fl_m2[idx], abs(fl_m2[idx]-fl_m1[idx]))
end

# ============================================================
#  2c. Glucose-maltose puls scenario (Model 2)
# ============================================================
println("\n=== 2c: Glucose-maltose puls scenario ===")

function run_glucose_maltose(N0, t_inf, moi, t_puls, malt_puls)
    p = make_p2(N0, t_inf, moi)
    u0 = zeros(23)   # Model 2 gebruikt 23-element vector
    u0[Model2.Sind]  = [4.44, 0.0, 5.42, 0.0]
    u0[Model2.Eind]  = [0.95, 0.01, 0.01, 0.01]
    u0[Model2.Nind]  = N0
    tspan = (0.0, p.duration)

    maltoseCondition(u,t,integrator)   = t == t_puls
    maltoseAffect!(integrator)         = integrator.u[Model2.Sind[2]] += malt_puls
    maltoseCallBack = DiscreteCallback(maltoseCondition, maltoseAffect!)

    infectionCondition(u,t,integrator) = t == p.infection_time
    infectionAffect!(integrator)       = integrator.u[Model2.Pfind] = p.infection_dose
    infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

    fbaUpdateTimepoints = collect(0:1/60:p.duration)
    fbaUpdateCondition(u,t,integrator) = t in fbaUpdateTimepoints
    fbaAffect!(integrator)             = Model2.fbaUpdate!(integrator.u, p)
    fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

    domainCondition(u,t,integrator) = any(x->x<0.0, u)
    domainAffect!(integrator)       = Model2.enforcePositiveDomain!(integrator.u)
    domainCallBack = DiscreteCallback(domainCondition, domainAffect!)

    problem = DDEProblem(Model2.simulate_dFBA!, u0, (p,t)->u0, tspan, p)
    return solve(problem, MethodOfSteps(Tsit5()),
        verbose=false, reltol=1e-4, abstol=1e-6,
        tstops=sort(unique([t_puls; p.infection_time; fbaUpdateTimepoints])),
        callback=CallbackSet(domainCallBack, maltoseCallBack, infectionCallBack, fbaCallBack))
end

t_inf_glc      = 5.0
malt_puls      = 2.337
delta_t_values = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
max_lys_glc    = Float64[]

println("  Maltosepuls sweep (t_inf=$t_inf_glc h, MOI=2.0)")
for dt in delta_t_values
    t_puls = t_inf_glc - dt
    t_puls < 0.1 && continue
    sol_gm = run_glucose_maltose(1e9, t_inf_glc, 2.0, t_puls, malt_puls)
    if SciMLBase.successful_retcode(sol_gm)
        Lys_ts = [sol_gm.u[i][Model2.lind] for i in eachindex(sol_gm.u)]
        push!(max_lys_glc, maximum(Lys_ts))
        println("  dt=$dt h voor infectie: max lysogenen=$(round(max_lys_glc[end], sigdigits=2))")
    else
        push!(max_lys_glc, 0.0)
        println("  dt=$dt h: simulatie gefaald")
    end
end

sol_std = Model2.run(make_p2(1e9, t_inf_glc, 2.0))
Lys_std = SciMLBase.successful_retcode(sol_std) ?
    maximum([sol_std.u[i][Model2.lind] for i in eachindex(sol_std.u)]) : 0.0

fig2c_left = bar(string.(delta_t_values[1:length(max_lys_glc)]), max_lys_glc,
    title="2c: Max lysogene cellen vs puls-timing",
    xlabel="Δt maltosepuls voor infectie [h]", ylabel="Max lysogene cellen/L",
    label="Glc+Malt puls", color=:teal, legend=:topright)
hline!(fig2c_left, [Lys_std], label="Standaard (altijd maltose)",
    color=:black, lw=2, linestyle=:dash)

idx_best = argmax(max_lys_glc)
dt_best  = delta_t_values[idx_best]
sol_best = run_glucose_maltose(1e9, t_inf_glc, 2.0, t_inf_glc-dt_best, malt_puls)
t_b   = sol_best.t
mal_b = [sol_best.u[i][Model2.Sind[2]] for i in eachindex(sol_best.u)]
N_b   = [sol_best.u[i][Model2.Nind]   for i in eachindex(sol_best.u)]
Lys_b = [sol_best.u[i][Model2.lind]   for i in eachindex(sol_best.u)]
Pf_b  = [max(sol_best.u[i][Model2.Pfind], 1.0) for i in eachindex(sol_best.u)]

fig2c_right = plot(t_b, [N_b Lys_b Pf_b],
    label=["Naïef N" "Lysogeen l" "Vrije fagen"],
    color=[:blue :green :red], lw=2, yscale=:log10, ylims=(1,:auto),
    title="2c: Tijdsreeks (puls $dt_best h voor infectie)",
    xlabel="t [h]", ylabel="cellen of fagen/L")
plot!(twinx(), t_b, mal_b, label="Maltose [mmol/L]",
    color=:orange, lw=2, linestyle=:dash,
    ylabel="Maltose [mmol/L]", legend=:topright)

fig2c = plot(fig2c_left, fig2c_right, layout=(1,2), size=(1000,420))
savefig(fig2c, "h2c_glucose_maltose_puls.png")
println("  Figuur opgeslagen: h2c_glucose_maltose_puls.png")

println("\n=== Hoofdstuk 2 voltooid ===")