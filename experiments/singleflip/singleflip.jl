# singleflip.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


function generate_single_flip_energy_runs(lattice::Lattice, copies::Int64, beta_values::Vector{Float64})
    tau_values = lattice.tau_values
    println(tau_values)
    energy_runs = [generate_energies(lattice, beta_values,1, tau_values, "single flip") for _ in 1:copies]
    return energy_runs
end

 
N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)
beta_values = 1 ./ generate_T_intervals(4.0, 0.25, 100)
copies = 20


datafile = "single_flips_lowtau_many_runs.csv"
folder = ""


time = now()
println("started")
results = generate_single_flip_energy_runs(lattice, copies, beta_values)
writedlm(joinpath(folder, datafile), results, ',')
println(now() - time)


