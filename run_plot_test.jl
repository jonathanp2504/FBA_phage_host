using Plots
println("PLOTS_OK")
# quick test plot
p = plot(rand(10))
try
    savefig(p, joinpath(@__DIR__, "test_plot.png"))
    println("SAVED test_plot.png")
catch e
    println("SAVE_FAILED: ", e)
end
