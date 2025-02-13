#saddles.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates
using Profile, ProfileView

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)
beta_values = 1 ./ generate_T_intervals(10.0, 0.25, 200)


datafile = "saddles.csv"
folder = "experiments\\saddles"
time = now()
results = generate_saddles(lattice, beta_values, 5, "k chain flip",20)
writedlm(joinpath(folder, datafile), results, ',')
println(now() - time)
p1 = plot()
scatter!(p1, results[2], results[1], markerstrokedwidth = 1, markersize = 1, label = "")
savefig(p1, joinpath(folder, "saddles.png"))

p2 = plot()
scatter!(p2, results[3], results[1], markerstrokedwidth = 1, markersize = 1, label = "")
savefig(p2, joinpath(folder, "saddles_beta.png"))