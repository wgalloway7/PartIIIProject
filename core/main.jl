# main.jl
using Random
using Plots

include("lattice.jl")
include("montecarlo.jl")

N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

function magnetisation(lattice::Lattice)
    return sum(lattice.grid) / (lattice.N * lattice.N)
end


function run_monte_carlo_for_beta(lattice::Lattice, beta::Float64, maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1), verbose::Bool=false)
    (converged, _, _, _) = run_metropolis_algorithm(lattice, beta, maximum_iterations, configuration_correlation_convergence_criteria, verbose=verbose)
    avg_magnetisation = magnetisation(lattice)
    return avg_magnetisation
end

function monte_carlo_vs_beta(lattice::Lattice, beta_range::Vector{Float64})
    avg_magnetisations = Float64[]
    
    for beta in beta_range
        avg_magnetisation = abs(run_monte_carlo_for_beta(lattice, beta))
        push!(avg_magnetisations, avg_magnetisation)
    end

    plot(beta_range, avg_magnetizations, xlabel="Beta", ylabel="Average Magnetization", label="Avg Magnetization", title="Average Magnetization vs Beta")
end

beta_range = collect(0.001:0.001:0.6)


monte_carlo_vs_beta(lattice, beta_range)
