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
k_values = [i for i in 1:9]
copies = 10
folder = "experiments\\k chain"

move = "k chain flip"

#n_multiplier = [2.0^(i-2) for i in 1:4]
n_multiplier = [1.0]
decorrelation_copies = 10
maximum_iterations = 1000

file_name_root =  "k chain flip"
data_file_root = "k chain flip"


for n in n_multiplier
    time_now = now()
    file_name = file_name_root * "n = $n.png"
    data_file = data_file_root * "n = $n.csv"

    figure_E_anneal(lattice, beta_values, k_values, copies, file_name, data_file, folder, N, move, decorrelation_copies,n,maximum_iterations)
    println(now() - time_now)
end



println("Finished all")
