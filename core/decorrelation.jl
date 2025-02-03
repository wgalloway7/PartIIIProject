#decorrelation.jl
using Random
using Plots
using Statistics
include("lattice.jl")
include("montecarlo.jl")

N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)

function generate_T_intervals(T_start::Float64, T_end::Float64, num_intervals::Int64)
    return T_start .* (T_end / T_start) .^ (range(0, 1, length=num_intervals))
end

function prepare_lattice!(lattice::Lattice,k::Int64 = 1; maximum_iterations::Int64=1000)
    lattice.grid = solved_configuration(lattice.N)
    run_metropolis_algorithm(lattice, 0.0, k, maximum_iterations = maximum_iterations)
end

function generate_decorrelation_n(lattice::Lattice, beta_values::Vector{Float64}; k::Int64 = 1, move::String = "single flip")
    prepare_lattice!(lattice, k)
    decorrelation_n = Int64[]
    for beta in beta_values
        (_, _, current_iteration, _) = run_metropolis_algorithm(lattice, beta, k, move, use_correlation = true, maximum_iterations = 10000)
        push!(decorrelation_n, current_iteration)
    end
    return decorrelation_n
end


function generate_energies(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, n_correlation::Vector{Int64}, move::String = "single flip")
    #n_correlation = generate_decorrelation_n(lattice, beta_values; k=k, move="single flip")
    energies = Float64[]
    for (i,beta) in enumerate(beta_values)
        run_metropolis_algorithm(lattice, beta, k, move, maximum_iterations = n_correlation[i], use_correlation = false)
        E = energy(lattice)
        push!(energies, E)
    end
    return energies
end

function concatenate_energies(lattice::Lattice, copies::Int64, beta_values::Vector{Float64}, k::Int64, move::String = "single flip")
    n_correlation = generate_decorrelation_n(lattice, beta_values; k=k, move=move)
    p2 = plot(beta_values, n_correlation)
    display(p2)
    energy_runs = [generate_energies(lattice, beta_values, k, n_correlation, move) for _ in 1:copies]
    return mean(hcat(energy_runs...), dims=2)[:]
end


beta_values = 1 ./ generate_T_intervals(0.01,5.0, 100)

#k_values = [i for i in 1:5]
k_values = [1]

copies = 1
filename = "test2.png"


#colors = cgrad(:RdBu, length(k_values))  # Set1 is a color palette with distinct colors
#p = plot()
#plot!(p, background_color = "#333333", gridcolor = :white, legend=:topright)
for (i, k) in enumerate(k_values)
    avg_energies = concatenate_energies(lattice, copies, beta_values, k, "single flip")
    #plot!(p, beta_values, avg_energies, label = "k = $k", color = colors[i])
end
    
#xlabel!("Beta", xlabelcolor = :white)
#ylabel!("Energy", ylabelcolor = :white)
#title!("Average Energy vs Beta for various k", titlecolor = :white)
    
# Set legend text to white
savefig(filename)
    
