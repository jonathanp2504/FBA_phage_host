using COBREXA, HiGHS, AbstractFBCModels, SBMLFBCModels
model_path = joinpath(@__DIR__, "iJO1366.xml")
@assert isfile(model_path)
m = convert(AbstractFBCModels.CanonicalModel.Model, load_model(model_path))
println("loaded OK")
