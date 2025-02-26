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
    plot!(p, beta_values, log10.(avg_correlation))
    
    xlabel!(p, "Beta")
    ylabel!(p, "log(Average n_correlation)")
    title!(p, "N = $N, $copies runs, $monte_carlo_steps  Monte Carlo steps")
    savefig(p, filename)
end


copies = 5
beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 100)
N_values = [10,20,50,100]
monte_carlo_steps = 10
lattice = Lattice(10)
for N in N_values
    time = now()
    println(N)
    lattice.N = N
    lattice.grid = solved_configuration(N)
    filename = joinpath("new_decorrelation_avg_N$N.png")
    figure_n_correlation(lattice, beta_values, copies, filename, N, monte_carlo_steps)
    println(now() - time)
    
end

monte_carlo_steps_values = [1,2,5,10,20]
lattice.N = 50
lattice.grid = solved_configuration(50)
for m in monte_carlo_steps_values
    time = now()
    println(m)
    filename = joinpath("new_decorrelation_avg_m$m.png")
    figure_n_correlation(lattice, beta_values, copies, filename, 50, m)
    println(now() - time)
end
    
    