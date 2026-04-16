using Pkg
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()
Pkg.precompile()
import Pkg
Pkg.add("DelayDiffEq")
Pkg.add("SciMLBase")
Pkg.add("JuMP")
using DelayDiffEq
using SciMLBase
using JuMP

