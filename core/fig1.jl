# main.jl
using Random
using Plots
using Statistics


include("lattice.jl")
include("montecarlo.jl")

N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

function run_monte_carlo_for_beta(lattice::Lattice, beta::Float64, k::Int64, copies::Int64 = 1; maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1), verbose::Bool=false, move::String = "single flip")
    energies = Float64[]
    for i in 1:copies
        (converged, _, _, _) = run_metropolis_algorithm(lattice, beta, k, move; maximum_iterations=maximum_iterations, configuration_correlation_convergence_criteria=configuration_correlation_convergence_criteria, verbose=verbose)
        E = energy(lattice)
        push!(energies, E)
    end
    
    return mean(energies)
end

function monte_carlo_vs_inverse_beta(lattice::Lattice, inverse_beta_range::Vector{Float64}, k::Int64, copies::Int64, move::String)
    avg_energies = Float64[]
    for inv_beta in inverse_beta_range
        beta = 1 / inv_beta
        avg_energy = abs(run_monte_carlo_for_beta(lattice, beta, k,copies; move=move))
        push!(avg_energies, avg_energy)
    end
    return avg_energies
    #plot(inverse_beta_range, avg_energies, xlabel="1/Beta", ylabel="Energy", label="Avg Energy")
end

#inverse_beta_range = collect(1.0:0.1:25)  # Adjust the range as needed
#k = 3

#plot(inverse_beta_range, monte_carlo_vs_beta(lattice, inverse_beta_range, k, 10, move = "k_chain_flip"), xlabel="1/Beta", ylabel="Energy", label="Avg Energy", title = "k chain")



function plot_inverse_beta_vs_energy(lattice::Lattice, inverse_beta_range::Vector{Float64}, k_values:: Vector{Int64}, copies::Int64, move::String, filename:: String)
    plot()
    for k in k_values
        avg_energies = monte_carlo_vs_inverse_beta(lattice, inverse_beta_range, k , copies, move)
        plot!(inverse_beta_range, avg_energies, label = "k = $k")
    end
    xlabel!("1/Beta")
    ylabel!("Energy")
    title!("Average Energy vs 1/Beta for various k")
    savefig(filename)
end

inverse_beta_range = collect(0.001:0.1:6)
k_values = [i for i in 1:15]
copies = 20
move = "k chain flip"
filename = "inverse_beta_vs_energy.png"

plot_inverse_beta_vs_energy(lattice, inverse_beta_range, k_values, copies, move, filename)