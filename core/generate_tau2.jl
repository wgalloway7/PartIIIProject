using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
include("../core/lattice.jl")
include("../core/montecarlo.jl")
include("../core/main.jl")


N = 20
res = 100
lattice = Lattice(N)
#beta_values= 1 ./ generate_T_intervals(10.0,0.2,res)
beta_values = 1 ./ collect(range(0.2,stop = 5.0, length=res))
time = now()
filename = "autocorrelation_short.csv"
tau_values = generate_tau_quick2(lattice, beta_values, true, filename = filename)
println(now() - time)
writedlm("tau2_values.csv", (beta_values, tau_values), ',')
