# singleflip.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


function generate_single_flip_energy_runs(lattice::Lattice, copies::Int64, beta_values::Vector{Float64}, monte_carlo_timesteps::Int64)
    n_correlation = generate_decorrelation_n(lattice, beta_values; k = 1, move = "single flip", maximum_iterations = lattice.N^2 * monte_carlo_timesteps, copies = 10)
    energy_runs = [generate_energies(lattice, beta_values, 1, n_correlation, "single flip") for _ in 1:copies]
    return energy_runs
end

 
N = 200
lattice = Lattice(N)
lattice.grid = solved_configuration(N)
beta_values = 1 ./ generate_T_intervals(10.0, 0.25, 100)
copies = 25
monte_carlo_timesteps = 20


datafile = "single_flips_N200.csv"
folder = ""



time = now()
println("started")
results = generate_single_flip_energy_runs(lattice, copies, beta_values, monte_carlo_timesteps)
writedlm(joinpath(folder, datafile), results, ',')
println(now() - time)


