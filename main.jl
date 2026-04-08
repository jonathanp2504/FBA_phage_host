using COBREXA
using HiGHS
using AbstractFBCModels
import SBMLFBCModels
using Pkg
#Pkg.add("OrdinaryDiffEq")
#Pkg.add("DifferentialEquations")
#Pkg.add("DelayDiffEq")
using Plots
using DifferentialEquations
using OrdinaryDiffEq

# 1. Model laden
model_path = joinpath(@__DIR__, "iJO1366.xml")
@assert isfile(model_path) "iJO1366.xml does not exist at path: $model_path"
model = convert(AbstractFBCModels.CanonicalModel.Model, load_model(model_path))


model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0 # bacterial growth is constrained!!!
model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].upper_bound = 1.0 # [1/h] moet bovenlimiet op staan anders te hoog
# 2. Cybernetische Parameters 
alpha_syn = 0.2*10   # Snelheid van enzymsynthese waardes van deze 2 nog opzoeken (stonden te laag)
beta_deg = 0.05*10   # Snelheid van enzymdegradatie
K_s = [0.0278, 0.0146, 0.0543, 0.0833]      # Affiniteit uit Tabel 7 (mmol/L)
tau = 1.0          # Latente periode (uren)
b = 170.0          # Burst size
alfa_ads = 1e-10  # Adsorptieconstante (L/gDW/h) lager dan bij Luan want eigenlijk meer deze waarde voor lambda bij te hoge alfa numerieke problemen), ik wil 2-step want moet voor hoge MOI zoals bij ons uiteindelijk gaat zijn
p_pref = [0.8925, 0.08925,  0.008925,  0.008925] # Voorkeurshiërarchie uit thesis
V_max = [15.0, 13.0, 11.0, 4.0]           # Opnamesnelheden (mmol/hr^-1)
E_coli_cellDW = 2.8e-13 # gDW per cel
infection_time = 5.0 # Tijdstip van infectie (uren)
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e", "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e", "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e", "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e", "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release = 5.66e-12 # glucose vrijgegeven bij lysis (mmol/burst) waarde gebaseerd op Luan table 7
duration = 15.0
n_hill = 2.0
# K_mu = 0.2
f_res = 0.05
K_mal = 0.01
mu_max_glc = 1.33
mu_max_mal = 1.26
mu_max_gly = 1.10
mu_max_ac = 0.29
# --- 2. De Struct aanmaken ---
# Zorg dat de volgorde in je Parameters-file exact matcht met deze aanroep:
include("./parameters.jl")

# Ladder configuratie
max_phages = 10 

# Definieer de indices voor de state vector u
# 1:4   -> Substraten
# 5:8   -> Enzymen
# 9:19  -> S_ladder (S0, S1, S2, ..., S10) -> 11 variabelen
# 20    -> I (Geïnfecteerd)
# 21    -> L (Lysogeen)
# 22    -> P (Vrije fagen)
# 23    -> Benz (Benzonase)

p = Parameters(
    # Cybernetica
    alpha_syn, 
    beta_deg, 
    K_s, 
    V_max, 
    p_pref, 
    
    # Phage-Host
    tau, 
    b, 
    E_coli_cellDW, 
    MW_values, 
    h_release, 

    # Model IDs
    "R_BIOMASS_Ec_iJO1366_core_53p95M", 
    exchange_ids, 
    all_ex_ids, 
    essentials_ids, 
    model, 
    0.0,            # mu start op 0
    zeros(4),       # q start op 0

    # Adsorption Ladder (Sequential Model)
    max_phages,     # 10
    1e-10,          # k_on (jouw alfa_ads)
    10.0,           # k_off (1/h) -> Faag valt relatief snel los
    5.0,            # k_ins (1/h) -> Gemiddeld 12 min voor injectie
    K_mal,          # 0.01 (Verzadiging van LamB)

    # Indices
    1:4,            # ind_subs
    5:8,            # ind_e
    9:19,           # ind_S_ladder
    20,             # ind_I
    21,             # ind_L
    22,             # ind_P
    infection_time, 

    # Benzonase
    23,             # ind_Benz
    "R_BENZ_prod", 
    0.05,           # k_tox
    0.1,            # beta_benz
    0.0,            # q_benz start op 0    
    
    # Luan Referenties
    mu_max_glc, 
    mu_max_mal, 
    mu_max_gly, 
    mu_max_ac)

include("./dFBA_function.jl")
# 5. INITIALISATIE
# [Glc, Mal, Glyc, Ac, e_glc, e_mal, e_Glyc, e_ac, S, I, L, P, Benz] (subs in mmol/l)
u0 = zeros(23)
u0[1:4] = [4.44, 2.337, 5.42, 0.0]    # Subs
u0[5:8] = [0.95, 0.01, 0.01, 0.01]    # Enzymen
u0[9]   = 1e6                         # S0 (alle cellen beginnen zonder fagen)
u0[10:19] .= 0                         # S1 t/m S10
u0[20:23]  .= 0                         # I
# Let op: check in je iJO1366.xml of het "atp_c" of "atp[c]" is!
benz_stoich = Dict(
    "M_atp_c" => 0, 
    "M_h2o_c" => 0, 
    "M_adp_c" =>  0, 
    "M_pi_c"  =>  0, 
    "M_h_c"   =>  0
)

# 2. Bouw het officiële Reaction object
benz_reaction = AbstractFBCModels.CanonicalModel.Reaction(
    stoichiometry = benz_stoich,
    lower_bound = 0.0,   # Wordt door fbaUpdate aangepast
    upper_bound = 1000.0
)

# 3. Voeg het toe aan het model dictionary
model.reactions["R_BENZ_prod"] = benz_reaction
fbaUpdate!(u0, p)
tspan = (0.0, duration) 

infectionCondition(u, t, integrator) = t == p.infection_time #pathway in FBA bouwen die iets zegt over hoeveel product er wordt gemaakt --> dan uw model voor set van verschillende infection times runnen en zien wat maximale productie geeft
infectionAffect!(integrator) = integrator.u[p.ind_P] = 1000.0
infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

fbaUpdateTimepoints = collect(0:0.1:10.0) # Elke 0.1 uur een FBA-update uitvoeren wordt gestored in fbaModel in parameters (sneller dan als we het niet zo doen) --> kleine bug glucose gaat onder 0 dus zet constraint op 0
fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
fbaAffect!(integrator) = fbaUpdate!(integrator.u, p)
fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

prob = DDEProblem(dFBA_phage_system, u0, (p,t)->u0, tspan, p)
sol = solve(prob, MethodOfSteps(Tsit5()), 
            reltol=1e-3, # Zet deze weer iets strenger voor detail
            abstol=1e-5,
            maxiters=1e5,
            tstops=[p.infection_time; fbaUpdateTimepoints],
            callback=CallbackSet(infectionCallBack, fbaCallBack),
            isoutofdomain=(u,p,t) -> any(x -> x < -1e-8, u))
            #met isoutofdomain stopt solver telkens na 6 uur omdat glucose onder 0 gaat, dus nu met discrete callback die elke 0.1 uur FBA update doet en daarin zetten we constraint op glucose op 0 als die onder 0 gaat, werkt goed en sneller dan zonder discrete callback (omdat we niet telkens opnieuw moeten oplossen maar gewoon de constraints aanpassen in het model)   
        #sol = solve(prob, MethodOfSteps(Rosenbrock23(autodiff=false)), 
            #reltol=1e-3, 
            #abstol=1e-5,
            #maxiters=1e6, # Verhoogd om niet te vroeg te stoppen
            #tstops=[p.infection_time; fbaUpdateTimepoints],
            #callback=CallbackSet(infectionCallBack, fbaCallBack),
            #isoutofdomain=(u,p,t)->any(x->x<0, u)) # Voorkomt negatieve populaties 
# Gebruik een stijvere solver met betere stabiliteit
          
# 6. UITGEBREID PLOTTEN
# Plot A: Substraatverloop
p1 = plot(sol, idxs=p.ind_subs, title="Substraten (mmol/L)", 
          label=["Glucose" "Maltose" "Glycerol" "Acetaat"], lw=2, xlabel= "t [h]", ylabel="Conc.")

# Plot B: Enzymatische Adaptatie 
p2 = plot(sol, idxs=p.ind_e, title="Enzym-niveaus (Cybernetische e)", 
          label=["e_Glc" "e_Mal" "e_Glyc" "e_Ac"], lw=2, ls=:dash, xlabel= "t [h]", ylabel="Relatief niveau")


# Plot C: Host populatie
# Bereken data voor de plots
S_totaal_data = [sum(u[p.ind_S_ladder]) for u in sol.u]
X_totaal_data = [getTotalBiomass(u, p) for u in sol.u]

p3 = plot(sol.t, X_totaal_data, title="Bacteriële Populatie", yscale=:log10, ylims=(1,:auto),
          label="Totaal X", color=:black, lw=3, ylabel="cells [1/L]")

# Hier plotten we de som van S0 t/m S10
plot!(p3, sol.t, S_totaal_data, label="Totaal Vatbaar (S)", color=:blue, lw=2)
plot!(p3, sol, idxs=[p.ind_I, p.ind_L], label=["Infect (I)" "Lysogeen (L)"], alpha=0.7)

# Plot D: Faag Lambda Dynamiek
p4 = plot(sol.t, [getPhages(u, p) for u in sol.u], title="Fagen (P)", yscale=:log10, ylims=(1,:auto),
          color=:red, lw=2, label="Phage Lambda", xlabel= "t [h]", ylabel="phages [1/L]")


# Plot E: Benzonase Accumulatie
p5 = plot(sol.t, [u[p.ind_Benz] for u in sol.u], 
          title="Benzonase Productie", 
          color=:green, lw=2, 
          xlabel="t [h]", ylabel="Conc. [mmol/L]",
          label="Benzonase")

# Plot voor de individuele ladder-treden
p6 = plot(title="Adsorptie Ladder Verdeling", yscale=:log10, 
          xlabel="t [h]", ylabel="cells [1/L]", ylims=(1, :auto))

# We loopen door de treden om ze allemaal toe te voegen
for i in 0:p.max_phages
    # i=0 is S0, i=1 is S1, etc.
    label_name = "S$i"
    # Gebruik een kleurverloop van blauw naar rood voor de bezettingsgraad
    plot!(p6, sol, idxs=p.ind_S_ladder[i+1], label=label_name, alpha=0.6)
end

# Combineer alle plots in een nieuwe layout
plot(p1, p2, p3, p4, p5, p6, 
                  layout=(3, 2), 
                  size=(1000, 1200), # Iets groter gemaakt voor de leesbaarheid
                  margin=5Plots.mm)
#volgende stap is nu productie erin te krijgen en in uw FBA inbouwen door bijvoorbeeld zelf reacties toe te voegen zoals in stoich matrix (bijvoorbeeld met insuline/benzoase (enzym misschien goede eerste stap, moet niet per se aan FBA omdat enzym) --> welke enzymes etc zijn nodig in E. Coli hiervoor en welke substraten moeten voor deze productie worden opgenomen, zie paper in edge)
#branch maken in github ook dat je terugkan als het niet lukt
# als ik dat heb zoek op optimization.jl (zoek package en neem door)