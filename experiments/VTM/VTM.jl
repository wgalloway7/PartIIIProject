# VTM.jl
# run to generate graphs for decorrelation curves for various beta
# TO DO: add csv generation

using Random
using Plots
using Statistics
using DelimitedFiles

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


N = 50
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

beta_values = 1 ./ generate_T_intervals(10.0, 0.8, 6)
filename = "hmmm.png"
m = 10

figure_correlation_decay(lattice, beta_values, 1, filename, 10000, m, 10000, "single flip")