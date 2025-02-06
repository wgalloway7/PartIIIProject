#main.jl
using Random
using Plots
using Statistics
include("lattice.jl")
include("montecarlo.jl")

N = 10
lattice = Lattice(N)
lattice.grid = solved_configuration(N)



function concatenate_energies(lattice::Lattice, copies::Int64, beta_values::Vector{Float64}, k::Int64, move::String = "single flip")
    n_correlation = generate_decorrelation_n(lattice, beta_values; k=k, move=move, maximum_iterations = 100000)
    energy_runs = [generate_energies(lattice, beta_values, k, n_correlation, move) for _ in 1:copies]
    return mean(hcat(energy_runs...), dims=2)[:]
end



function figure_E_anneal(lattice::Lattice, beta_values::Vector{Float64}, k_values::Vector{Int64}, copies::Int64, filename::String, N::Int64, move::String = "single flip")
    colors = cgrad(:RdBu, length(k_values))  # Set1 is a color palette with distinct colors
    p = plot()
    plot!(p, background_color = "#333333", gridcolor = :white, legend=:topright)
    for (i, k) in enumerate(k_values)
        avg_energies = concatenate_energies(lattice, copies, beta_values, k, move)
        plot!(p, beta_values, avg_energies, label = "k = $k", color = colors[i])
    end

    xlabel!(p,"Beta", xlabelcolor = :white)
    ylabel!(p,"Energy", ylabelcolor = :white)
    title!(p,"N = $N, copies = $copies, $move", titlecolor = :white)
    savefig(p, filename)
end

beta_values = 1 ./ generate_T_intervals(10.0,0.5, 100)

#k_values = [i for i in 1:5]
k_values = [1]

copies = 1
filename = "test2.png"

figure_E_anneal(lattice, beta_values, k_values, copies, filename,N,"single flip")
