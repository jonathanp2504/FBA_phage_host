using COBREXA
import SBMLFBCModels
using HiGHS
using AbstractFBCModels
using ConstraintTrees

### Documentatie omtrent de FBA library!
# https://cobrexa.github.io/COBREXA.jl/stable/

if !isfile("e_coli_core.xml")
    download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.xml",
    "e_coli_core.xml",
    "b4db506aeed0e434c1f5f1fdd35feda0dfe5d82badcfda0e9d1342335ab31116",)
end



mutable struct Metabolism
    growthRate::Float64
    fluxes::ComponentArray{Float64}

    fbaModel::AbstractFBCModels.CanonicalModel.Model

    function Metabolism(modelFile::String)
        fbaModel = convert(AbstractFBCModels.CanonicalModel.Model, load_model(modelFile));
        fbaModel.reactions["R_BIOMASS_Ecoli_core_w_GAM"].upper_bound = 2.0 
        return new(0.0, ComponentArray{Float64}(R_EX_glc__D_e=0.0, R_EX_mal__L_e=0.0), fbaModel)
    end
end


function updateMetabolism!(metabolism::Metabolism, medium::ComponentArray{Float64})
    try
        for key in keys(u.medium)
            metabolism.fbaModel.reactions[String(key)].lower_bound = -medium[key]
        end
        fbaSol = flux_balance_analysis(metabolism.fbaModel, optimizer = HiGHS.Optimizer)
        for key in keys(medium)
            metabolism.fluxes[String(key)] = fbaSol.fluxes[String(key)]
        end
        metabolism.growthRate = fbaSol.fluxes["R_BIOMASS_Ecoli_core_w_GAM"]
    catch
        metabolism.growthRate = 0.0
        metabolism.fluxes .= 0.0
    end
end

function getGrowthRate(metabolism::Metabolism)
    return metabolism.growthRate
end

function getFluxes(metabolism::Metabolism)
    return metabolism.fluxes
end
