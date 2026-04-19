using Pkg
Pkg.activate(".")
Pkg.add([
    "OrdinaryDiffEq",
    "DelayDiffEq",
    "SciMLBase",
    "OrdinaryDiffEqCore",
    "UnPack",
    "Plots",
    "COBREXA",
    "AbstractFBCModels",
    "SBMLFBCModels",
    "JuMP",
    "HiGHS"
])
Pkg.instantiate()
Pkg.precompile()