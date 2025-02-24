using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

data = readdlm("VTM_fit_kchain98.csv",',')

tau = data[1,:]
beta = data[3,:]

println(sizeof(tau))
println(sizeof(beta))
println(typeof(tau))
println(typeof(beta))

VFT_fit(beta, tau, "VFTfit98.png", 100, true)