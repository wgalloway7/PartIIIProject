#old program
# can just do include main.jl
# generates decorrelation numbers, produces csv and figure
using Random
using Plots
using Statistics
include("lattice.jl")
include("montecarlo.jl")

N = 10
lattice = Lattice(N)

beta_values = 1 ./ generate_T_intervals(10.0, 0.05, 100)
k_values = [1]
m = 10  # Number of runs to average over

function figure_n_correlation(lattice::Lattice, beta_values::Vector{Float64}, k_values::Vector{Int64}, m::Int64, filename::String, N::Int64)
    colors = cgrad(:RdBu, length(k_values))
    p = plot()
    plot!(p, background_color="#333333", gridcolor=:white, legend=:topright)
    
    for (i, k) in enumerate(k_values)
        all_correlations = []
        
        for _ in 1:m
            push!(all_correlations, generate_decorrelation_n(lattice, beta_values; k=k, move="single flip", maximum_iterations=10000))
        end
        
        avg_correlation = mean(all_correlations)
        plot!(p, beta_values, avg_correlation, label="k = $k", color=colors[i])
    end
    
    xlabel!(p, "Beta", xlabelcolor=:white)
    ylabel!(p, "Average n_correlation", ylabelcolor=:white)
    title!(p, "Average n_correlation vs Beta for various k, N = $N, number of runs = $m", titlecolor=:white)
    savefig(p, filename)
end

filename = "decorrelation_avg"
figure_n_correlation(lattice, beta_values, k_values, m, filename, N)