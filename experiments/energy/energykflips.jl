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

N = 20
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 100)
k_values = [i for i in 1:7]
copies = 100
root = "kenergy"
folder = "experiments\\energy"


moves = ["k chain flip", "k line flip", "unconstrained k flip"]

for move in moves
    file_name =  move* "_" * root * ".png"
    data_file = move * "_" * root * ".csv"
    figure_E_anneal(lattice, beta_values, k_values, copies, file_name, data_file, folder, N, move)
    println("Finished $move")
end
println("Finished all")
