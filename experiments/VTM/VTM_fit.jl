# VTM.jl
# run to generate fits for decorrelation curves for various beta

using Random
using Plots
using Statistics
using DelimitedFiles

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


N = 50
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

beta_values = 1 ./ generate_T_intervals(10.0, 0.8, 6)

m = 10

root = "VTM_fit"
folder = "experiments\\VTM"
moves = ["single flip","k chain flip", "k line flip", "unconstrained k flip"]

for move in moves
    filename =  move* "_" * root * ".png"
    datafile = move * "_" * root * ".csv"
    fit_VTM(lattice, beta_values, 1,filename, folder, datafile, "single flip", 2500, 10000, 10, true)
    println(move)
end

println("done")