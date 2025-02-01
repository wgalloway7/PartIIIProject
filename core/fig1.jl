# main.jl
using Random
using Plots
using Statistics


include("lattice.jl")
include("montecarlo.jl")

N = 50
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

function run_monte_carlo_for_beta(lattice::Lattice, beta::Float64, k::Int64, copies::Int64 = 1; maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1), verbose::Bool=false, move::String = "single flip")
    energies = Float64[]
    for i in 1:copies
        (converged, _, _, _) = run_metropolis_algorithm(lattice, beta, k, move; maximum_iterations=maximum_iterations, configuration_correlation_convergence_criteria=configuration_correlation_convergence_criteria, verbose=verbose)

        E = energy(lattice)
        push!(energies, E)
        lattice.grid = solved_configuration(N)
    end
    
    return mean(energies)
end

function monte_carlo_vs_beta(lattice::Lattice, beta_range::Vector{Float64}, k::Int64, copies::Int64, move::String)
    avg_energies = Float64[]
    for beta in beta_range

        avg_energy = abs(run_monte_carlo_for_beta(lattice, beta, k,copies; move=move))
        push!(avg_energies, avg_energy)
    end
    return avg_energies
end



function plot_beta_vs_energy(lattice::Lattice, beta_range::Vector{Float64}, k_values::Vector{Int64}, copies::Int64, move::String, filename::String)
    # Use a color gradient that ensures distinct colors
    colors = cgrad(:RdBu, length(k_values))  # Set1 is a color palette with distinct colors
    p = plot()
    plot!(p, background_color = "#333333", gridcolor = :white, legend=:topright)
    for (i, k) in enumerate(k_values)
        avg_energies = monte_carlo_vs_beta(lattice, beta_range, k , copies, move)
        plot!(p, beta_range, avg_energies, label = "k = $k", color = colors[i])
    end
    
    xlabel!("Beta", xlabelcolor = :white)
    ylabel!("Energy", ylabelcolor = :white)
    title!("Average Energy vs Beta for various k", titlecolor = :white)
    
    # Set legend text to white
    savefig(filename)
end

beta_range = collect(range(0.1, 3.0, length = 10))
k_values = [i for i in 1:20]
copies = 20
move = "unconstrained k flip"
filename = "beta_vs_energy.png"

plot_beta_vs_energy(lattice, beta_range, k_values, copies, move, filename)