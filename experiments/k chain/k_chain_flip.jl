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

N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 100)
k_values = [i for i in 1:49]
copies = 10
folder = ""

move = "k chain flip"

time = now()
file_name =  "bigger k chain flip.png"
data_file = "bigger k chain flip.csv"

decorrelation_copies = 10
time = now()
println("started")
figure_E_anneal(lattice, beta_values, k_values, copies, file_name, data_file, folder, N, move)
println(time - now())
println("Finished all")
