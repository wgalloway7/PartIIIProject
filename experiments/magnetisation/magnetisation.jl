using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
using Base.Threads
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)
beta_values = 1 ./ generate_T_intervals(4.0, 0.25, 100)
copies = 5

magnetisation_runs = []
move = "single flip"
for i in 1:copies
    println("Starting copy $(i)")
    time = now()
    mag_run = generate_magnetisation(lattice, beta_values, 1, lattice.tau_values, move)
    println("Finished copy $(i) in $(now() - time)")
    push!(magnetisation_runs, mag_run)
end
writedlm("magnetisation.csv", hcat(beta_values, magnetisation_runs...), ',')

p = plot()
for i in 1:copies
    plot!(p, beta_values, abs(magnetisation_runs[i]), label="Copy $i")
end
savefig(p, "magnetisation.png")

