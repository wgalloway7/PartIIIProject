# energyflips.jl
# run to generate csv and graphs for energy/beta for various k and different moves
# adjust copies (number of energy measurements),beta_values (temperature range) and lattice size N as needed
using Random
using Plots
using Statistics
using DelimitedFiles

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 10)
k_values = [i for i in 1:3]
copies = 1
root = "kenergy"
folder = "experiments\\energy"


moves = ["k chain flip", "k line flip", "unconstrained k flip"]
figure_E_anneal(lattice, beta_values, k_values, copies, filename, datafile, folder, N, "single flip")

for move in moves
    filename =  move* "_" * root * ".png"
    datafile = move * "_" * root * ".csv"
    figure_E_anneal(lattice, beta_values, k_values, copies, filename, datafile, folder, N, move)
end
