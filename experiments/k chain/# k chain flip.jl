# k chain flip
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
copies = 10
folder = "experiments\\k chain"

move = "k chain flip"


file_name =  "k chain flip.png"
data_file = "k chain flip.csv"
n_multiplier = 1

decorrelation_copies = 10
time = now()
figure_E_anneal(lattice, beta_values, k_values, copies, file_name, data_file, folder, N, move, decorrelation_copies, n_multiplier)
println(now() - time)
println("Finished all")
