# VTM.jl
# run to generate graphs for decorrelation curves for various beta
# TO DO: add csv generation

using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

beta_values = 1 ./ generate_T_intervals(10.0, 0.8, 10)
filename = "VTM_kchain22.png"
datafile = "VTM_kchain22.csv"
copies = 50

measurement_MC_steps = 10
cooling_MC_steps = 10
k = 22

println("started")
time = now()
figure_correlation_decay(lattice, beta_values, 51, filename, datafile, measurement_MC_steps, copies, cooling_MC_steps, "k chain flip")
println(now() - time)
println("done")