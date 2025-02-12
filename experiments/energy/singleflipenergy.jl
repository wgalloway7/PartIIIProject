# singleflipenergy.jl
# run to generate csv and graphs for single flips
# adjust copies (number of energy measurements),beta_values (temperature range) and lattice size N as needed
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.25, 10)
k_values = [1]
copies = 2
filename = "single_flip.png"
datafile = "single_flips.csv"
folder = "experiments\\energy"

time = now()
figure_E_anneal(lattice, beta_values, k_values, copies, filename, datafile, folder, N, "single flip")
println(now() - time)
