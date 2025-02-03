#main.jl
using Random
using Plots
using Statistics

include("lattice.jl")
include("montecarlo.jl")

N = 20
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

function generate_T_intervals(T_start::Float64, T_end::Float64, num_intervals::Int64)
    return T_start .* (T_end / T_start) .^ (range(0, 1, length=num_intervals))
end



function prepare_lattice!(lattice::Lattice,k::Int64 = 1; maximum_iterations::Int64=1000)
    lattice.grid = solved_configuration(lattice.N)
    run_metropolis_algorithm(lattice, 0.0, k, maximum_iterations = maximum_iterations)
end

function generate_energies(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, move::String = "single flip")
    # prepare the lattice
    prepare_lattice!(lattice, k)
    energies = Float64[]
    for beta in beta_values
        (converged, _, _, _) = run_metropolis_algorithm(lattice, beta, k, move)
        E = energy(lattice)
        push!(energies, E)
    end
    return energies
end

function concatenate_energies(lattice::Lattice, copies::Int64, beta_values::Vector{Float64}, k::Int64, move::String = "single flip")
    energy_runs = [generate_energies(lattice, beta_values, k) for _ in 1:copies]
    return mean(hcat(energy_runs...), dims=2)[:]
end

beta_values = 1 ./ generate_T_intervals(5.0,0.01, 50)
scatter(1 ./ beta_values)
scatter(beta_values)


k_values = [i for i in 1:5]
copies = 50
filename = "test.png"


colors = cgrad(:RdBu, length(k_values))  # Set1 is a color palette with distinct colors
p = plot()
plot!(p, background_color = "#333333", gridcolor = :white, legend=:topright)
for (i, k) in enumerate(k_values)
    avg_energies = concatenate_energies(lattice, copies, beta_values, k, "k chain flip")
    plot!(p, beta_values, avg_energies, label = "k = $k", color = colors[i])
end
    
xlabel!("Beta", xlabelcolor = :white)
ylabel!("Energy", ylabelcolor = :white)
title!("Average Energy vs Beta for various k", titlecolor = :white)
    
# Set legend text to white
savefig(filename)
    
    

