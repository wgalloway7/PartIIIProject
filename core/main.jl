using Random
using Plots
using Statistics
using DelimitedFiles
include("lattice.jl")
include("montecarlo.jl")

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

function concatenate_energies(lattice::Lattice, copies::Int64, beta_values::Vector{Float64}, k::Int64, move::String = "single flip")
    n_correlation = generate_decorrelation_n(lattice, beta_values; k=k, move=move, maximum_iterations = 100000)
    energy_runs = [generate_energies(lattice, beta_values, k, n_correlation, move) for _ in 1:copies]
    return mean(hcat(energy_runs...), dims=2)[:]
end

function figure_E_anneal(lattice::Lattice, beta_values::Vector{Float64}, k_values::Vector{Int64}, copies::Int64, filename::String, datafile::String, folder::String, N::Int64, move::String = "single flip")
    colors = cgrad(:RdBu, length(k_values))  # Set1 is a color palette with distinct colors
    p = plot()
    plot!(p, background_color = "#333333", gridcolor = :white, legend=:topright)
    
    all_data = []  # Store data for writing to file
    for (i, k) in enumerate(k_values)
        avg_energies = concatenate_energies(lattice, copies, beta_values, k, move)
        push!(all_data, avg_energies)
        plot!(p, beta_values, avg_energies, label = "k = $k", color = colors[i])
    end
    
    xlabel!(p,"Beta", xlabelcolor = :white)
    ylabel!(p,"Energy", ylabelcolor = :white)
    title!(p,"N = $N, copies = $copies, $move", titlecolor = :white)
    savefig(joinpath(folder, filename))
    
    # Save data to file
    writedlm(joinpath(folder, datafile), hcat(beta_values, all_data...), ',')
end


N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 10)
k_values = [i for i in 1:2]
copies = 2
filename = "test3.png"
datafile = "energies.csv"
folder = "energies"

figure_E_anneal(lattice, beta_values, k_values, copies, filename, datafile, folder, N, "single flip")
println("done")