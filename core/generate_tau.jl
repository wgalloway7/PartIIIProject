using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
include("../core/lattice.jl")
include("../core/montecarlo.jl")
include("../core/main.jl")


N = 100
res = 50
lattice = Lattice(N)
beta_values= 1 ./ generate_T_intervals(10.0,0.2,res)
time = now()
filename = "autocorrelation.csv"
tau_values = generate_tau_quick(lattice, beta_values, copies = 10, save_file = true, tau0_Tmin = 1000)
println(now() - time)
writedlm("tau_values.csv", (beta_values, tau_values), ',')
