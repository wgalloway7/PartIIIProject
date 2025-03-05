#old program
# can just do include main.jl
# generates decorrelation numbers, produces csv and figure
using Random
using Plots
using Statistics
using Dates
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")



function figure_n_correlation(lattice::Lattice, beta_values::Vector{Float64}, copies::Int64, filename::String, N::Int64, monte_carlo_steps::Int64)
    # currently takes the average of m runs
    # but generate_decorrelation_n now does that internally
    p = plot()
    #plot!(p, background_color="#333333", gridcolor=:white, legend=:topright)
    avg_correlation = generate_decorrelation_n(lattice, beta_values; k=1, move="single flip", maximum_iterations=monte_carlo_steps * N^2, copies = copies)
    plot!(p, beta_values, (avg_correlation / N^2))
    
    xlabel!(p, "Beta")
    ylabel!(p, "MC steps")
    title!(p, "N = $N, $copies runs, $monte_carlo_steps  Monte Carlo steps")
    savefig(p, filename)
end


copies = 1
beta_values = 1 ./ generate_T_intervals(10.0, 0.1, 2)
monte_carlo_steps = 1
N = 10
lattice = Lattice(N)

time = now()
println("started")
figure_n_correlation(lattice, beta_values, copies, "test.png", N, monte_carlo_steps)
println(now() - time)

# 100 beta values
# 5 copies
# 10000 monte carlo steps
# 100 lattice size
# 8935923 milliseconds (2 hours, 28 minutes, 19 seconds)
    