using Pkg
#Pkg.add("OrdinaryDiffEq")
#Pkg.add("DifferentialEquations")
#Pkg.add("DelayDiffEq")
using Plots
using COBREXA
using AbstractFBCModels
using DifferentialEquations
import COBREXA: StandardModel, Metabolite, Reaction
include("./parameters.jl")
include("./FBA.jl")
include("./dFBA_function.jl")

# 1. Model laden
model_path = joinpath(@__DIR__, "iJO1366.xml")
@assert isfile(model_path) "iJO1366.xml does not exist at path: $model_path"
# Gebruik deze dictionary voor je add_reaction!
benz_stoich = Dict(
    # Aminozuren (Inputs vanuit het Cytosol, sequentie via uniprot)
    "M_ala__L_c" => -33.0, 
    "M_arg__L_c" => -13.0, 
    "M_asn__L_c" => -22.0, 
    "M_asp__L_c" => -16.0, 
    "M_cys__L_c" => -4.0,  
    "M_gln__L_c" => -12.0, 
    "M_glu__L_c" => -10.0, 
    "M_gly_c"    => -21.0, 
    "M_his__L_c" => -4.0, 
    "M_ile__L_c" => -8.0,  
    "M_leu__L_c" => -21.0, 
    "M_lys__L_c" => -14.0, 
    "M_met__L_c" => -3.0,  
    "M_phe__L_c" => -8.0,  
    "M_pro__L_c" => -10.0, 
    "M_ser__L_c" => -19.0, 
    "M_thr__L_c" => -15.0, 
    "M_trp__L_c" => -5.0,  
    "M_tyr__L_c" => -10.0, 
    "M_val__L_c" => -10.0,

    # Energieverbruik volgens de paper (Totaal 266 AA)
    # 1. ATP Gedeelte (2 ATP per AA = 532)
    "M_atp_c"    => -532.0,
    "M_adp_c"    =>  532.0, # (Let op: Strikt genomen AMP, maar in FBA vaak ADP voor balans)
    
    # 2. GTP Gedeelte (2 GTP per AA = 532)
    "M_gtp_c"    => -532.0,
    "M_gdp_c"    =>  532.0,

    # Overige bijproducten voor de massa-balans (P_i en H+)
    "M_h2o_c"    => -1064.0,
    "M_pi_c"     =>  1064.0,
    "M_h_c"      =>  1064.0,
    # Product (De 'M_benzonase_c' die je zelf aanmaakt)
    "M_benzonase_c" => 1.0 
)
using AbstractFBCModels
import SBMLFBCModels

# 1. Laad het model (dit geeft een SBMLFBCModel terug)
raw_model = load_model(model_path)

# 2. Converteer het naar een CanonicalModel om het bewerkbaar te maken
# Dit lost de MethodError op
model = convert(AbstractFBCModels.CanonicalModel.Model, raw_model)

# 3. Maak je 'lege' objecten aan voor de toevoegingen
# We voegen ze direct toe aan de dictionaries van het geconverteerde model
model.metabolites["M_benzonase_c"] = AbstractFBCModels.CanonicalModel.Metabolite()
model.metabolites["M_benzonase_c"].name = "Benzonase"
model.metabolites["M_benzonase_c"].compartment = "c"

model.metabolites["M_benzonase_e"] = AbstractFBCModels.CanonicalModel.Metabolite()
model.metabolites["M_benzonase_e"].name = "Benzonase (extracellular)"
model.metabolites["M_benzonase_e"].compartment = "e"

model.reactions["R_BENZ_prod"] = AbstractFBCModels.CanonicalModel.Reaction()
model.reactions["R_BENZ_prod"].name = "Benzonase production"
model.reactions["R_BENZ_prod"].lower_bound = 0.0
model.reactions["R_BENZ_prod"].upper_bound = 1000.0
model.reactions["R_BENZ_prod"].stoichiometry = benz_stoich # Je Dict

# --- STAP C: Benzonase Export (NIEUW: van cel naar medium) ---
model.reactions["R_BENZ_export"] = AbstractFBCModels.CanonicalModel.Reaction()
model.reactions["R_BENZ_export"].name = "Benzonase secretion"
model.reactions["R_BENZ_export"].lower_bound = 0.0
model.reactions["R_BENZ_export"].upper_bound = 1000.0
model.reactions["R_BENZ_export"].stoichiometry = Dict(
    "M_benzonase_c" => -1.0, 
    "M_benzonase_e" => 1.0
)
# STAP D: De Exchange reactie
model.reactions["R_EX_benz_e"] = AbstractFBCModels.CanonicalModel.Reaction()
model.reactions["R_EX_benz_e"].name = "Benzonase exchange (Sink)"
model.reactions["R_EX_benz_e"].lower_bound = 0.0      # De cel kan het niet 'opeten' van buitenaf
model.reactions["R_EX_benz_e"].upper_bound = 1000.0   # Het kan onbeperkt wegstromen
model.reactions["R_EX_benz_e"].stoichiometry = Dict(
    "M_benzonase_e" => -1.0                           # Het metaboliet 'verdwijnt' hier uit de berekening
)
model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].lower_bound = 0.0 # bacterial growth is constrained!!!
#model.reactions["R_BIOMASS_Ec_iJO1366_core_53p95M"].upper_bound = 2.0 # [1/h] moet bovenlimiet op staan anders te hoog
# 2. Cybernetische Parameters 
alpha_syn = 0.2*10   # Snelheid van enzymsynthese waardes van deze 2 nog opzoeken (stonden te laag)
beta_deg = 0.05*10   # Snelheid van enzymdegradatie
K_s = [0.0278, 0.0146, 0.0543, 0.0833]      # Affiniteit uit Tabel 7 (mmol/L)
tau = 1.0          # Latente periode (uren)
b = 170.0          # Burst size
alfa_ads = 1e-10  # Adsorptieconstante (L/gDW/h) lager dan bij Luan want eigenlijk meer deze waarde voor lambda bij te hoge alfa numerieke problemen), ik wil 2-step want moet voor hoge MOI zoals bij ons uiteindelijk gaat zijn
p_pref = [0.8925, 0.08925,  0.008925,  0.008925] # Voorkeurshiërarchie uit thesis
V_max = [12.7, 3.75, 0.0, 4.0]           # Opnamesnelheden (mmol/hr^-1)
E_coli_cellDW = 1.0e-12 # gDW per cel
infection_time = 2.0 # Tijdstip van infectie (uren)
essentials_ids = ["R_EX_o2_e", "R_EX_nh4_e", "R_EX_pi_e", "R_EX_so4_e", "R_EX_k_e", "R_EX_mg2_e", "R_EX_ca2_e", "R_EX_cl_e", "R_EX_fe2_e", "R_EX_fe3_e", "R_EX_mn2_e", "R_EX_zn2_e", "R_EX_cu2_e", "R_EX_cobalt2_e", "R_EX_mobd_e", "R_EX_thi_e", "R_EX_ni2_e", "R_EX_sel_e", "R_EX_slnt_e", "R_EX_tungs_e"]
exchange_ids   = ["R_EX_glc__D_e", "R_EX_malt_e", "R_EX_glyc_e", "R_EX_ac_e"]
all_ex_ids     = [id for id in keys(model.reactions) if startswith(id, "R_EX_")]
MW_values      = [180.16, 342.3, 92.09, 60.05]
h_release = 0.0 #6.0e-12 # glucose vrijgegeven bij lysis (mmol/burst) waarde gebaseerd op Luan table 7
duration = 15.0
n_hill = 2.0
# K_mu = 0.2
f_res = 0.05
K_mal = 0.01
mu_max_glc = 1.33
mu_max_mal = 1.26
mu_max_gly = 1.10
mu_max_ac = 0.29
mu_max_vector = [1.33, 1.26, 1.10, 0.29] 
# 2. Bereken e_max vector (Steady-state: synthese / (degradatie + groei))
# Formule: (alpha_syn + delta) / (beta_deg + mu_max)
e_max_vector = (alpha_syn .+ 0.001) ./ (beta_deg .+ mu_max_vector)# glc, mal, glyc, ac
# --- 2. De Struct aanmaken ---
# Zorg dat de volgorde in je Parameters-file exact matcht met deze aanroep:
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
    0.0,               # mu_N (startwaarde groei naïef)
    zeros(4),          # q_N  (startwaarde opname naïef)
    0.0,               # mu_l (startwaarde groei lysogeen)
    zeros(4),          # q_l  (startwaarde opname lysogeen)
    0.0,               # q_benz_l (startwaarde productie benzonase)

    # Adsorption Ladder (Sequential Model)
    1e-10,          # k_on (jouw alfa_ads)
    10.0,           # k_off (1/h) -> Faag valt relatief snel los
    5.0,            # k_ins (1/h) -> Gemiddeld 12 min voor injectie
    K_mal,          # 0.01 (Verzadiging van LamB)

    infection_time, 

    
    "R_BENZ_prod", 
    0.05,           # k_tox
    0.1,            # beta_benz
    0.0,            # q_benz start op 0    
    mu_max_vector,  # DE NIEUWE VECTOR
    e_max_vector,    # DE NIEUWE VECTOR
    0.0015            #f_prod (% van de totale import wordt opgeofferd aan benzonase productie)
)

include("./dFBA_function.jl")
# 5. INITIALISATIE
# [Glc, Mal, Glyc, Ac, e_glc, e_mal, e_Glyc, e_ac, S, I, L, P, Benz] (subs in mmol/l)
u0 = zeros(23)
u0[Sind] = [4.44, 2.337, 5.42, 0.0]    # Subs
u0[Eind] = [0.95, 0.01, 0.01, 0.01]    # Enzymen
u0[Nind]   = 1e9                         # S0 (alle cellen beginnen zonder fagen)

fbaUpdate!(u0, p)
tspan = (0.0, duration) 

infectionCondition(u, t, integrator) = t == p.infection_time #pathway in FBA bouwen die iets zegt over hoeveel product er wordt gemaakt --> dan uw model voor set van verschillende infection times runnen en zien wat maximale productie geeft
infectionAffect!(integrator) = integrator.u[Pfind] = 1e6
infectionCallBack = DiscreteCallback(infectionCondition, infectionAffect!)

fbaUpdateTimepoints = collect(0:1/60:10.0) # Elke 0.1 uur een FBA-update uitvoeren wordt gestored in fbaModel in parameters (sneller dan als we het niet zo doen) --> kleine bug glucose gaat onder 0 dus zet constraint op 0
fbaUpdateCondition(u, t, integrator) = t in fbaUpdateTimepoints
fbaAffect!(integrator) = fbaUpdate!(integrator.u, p)
fbaCallBack = DiscreteCallback(fbaUpdateCondition, fbaAffect!)

# positive domain
domainCondition(u, t, integrator) = any(x -> x < 0.0, u)
domainAffect!(integrator) = map!(v -> max(v, 0.0), integrator.u)
domainCallBack = ContinuousCallback(domainCondition, domainAffect!)

prob = DDEProblem(dFBA_phage_system, u0, (p,t)->u0, tspan, p)
sol = solve(prob, MethodOfSteps(Tsit5()), 
            reltol=1e-4, # Zet deze weer iets strenger voor detail
            abstol=1e-6,
            maxiters=1e5,
            tstops=[p.infection_time; fbaUpdateTimepoints],
            callback=CallbackSet(infectionCallBack, fbaCallBack, domainCallBack))
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
include("./plotting.jl")
plotAll(sol, p)