using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
using Base.Threads

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

N =200
lattice = Lattice(N)
steps = 10
move = "single flip"

lattice.grid = random_configuration(N, 0.0)
p = plot()
heatmap!(p, lattice.grid)
savefig(p, "initial.png")

for i in 1:steps
    for i in 1:(N^2 / 5)
        monte_carlo_timestep!(lattice, move, 1.0)
    end
    heatmap!(p, lattice.grid)
    savefig(p, "step_$(i).png")
    
    
end