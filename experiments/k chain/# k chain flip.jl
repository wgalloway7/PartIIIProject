# k chain flip
# run to generate csv and graphs for energy/beta for various k and different moves
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


beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 10)
k_values = [1]
copies = 1
folder = "experiments\\k chain"

move = "single flip"


file_name =  "k chain flip.png"
data_file = "k chain flip.csv"
n_multiplier = 1
maximum_iterations = 100000

decorrelation_copies = 1
time = now()
figure_E_anneal(lattice, beta_values, k_values, copies, file_name, data_file, folder, N, move, decorrelation_copies, n_multiplier, maximum_iterations)
println(now() - time)
println("Finished all")
