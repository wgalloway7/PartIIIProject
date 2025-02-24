# k chain tau
# run to geniperate csv and graphs for energy/beta for various k and different moves
# adjust copies (number of energy measurements),beta_values (temperature range) and lattice size N as needed
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates


include("../../../core/lattice.jl")
include("../../../core/montecarlo.jl")
include("../../../core/main.jl")

function figure_E_vary_tau(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, copies::Int64, filename::String, datafile::String, folder::String, N::Int64, decorrelation_n_multiplier_values::Vector{Float64}, move::String = "single flip", decorrelation_copies::Int64 = 1, monte_carlo_timesteps::Int64 = 10)
    colors = cgrad(:RdBu, length(decorrelation_n_multiplier_values))  # Set1 is a color palette with distinct colors
    p = plot()
    plot!(p, background_color = "#333333", gridcolor = :white, legend=:topright)
    
    all_data = []  # Store data for writing to file
    #generate decorrelation n for single flip ie k =1
    n_correlation = generate_decorrelation_n(lattice, beta_values; k=1, move="single flip", maximum_iterations=monte_carlo_timesteps * lattice.N^2, copies =decorrelation_copies)
    println(n_correlation)
    scaled_n_correlation = Int64.(ceil.(n_correlation .* first(decorrelation_n_multiplier_values)))
    for (i, tau) in enumerate(decorrelation_n_multiplier_values)
        scaled_n_correlation = Int64.(ceil.(n_correlation .* tau))
        avg_energies = concatenate_energies(lattice, copies, beta_values, k, scaled_n_correlation, move)
        push!(all_data, avg_energies)
        plot!(p, beta_values, avg_energies, label = "N = $tau N_tau", color = colors[i])
    end
    avg_energies = concatenate_energies(lattice, copies, beta_values, 1, n_correlation, "single flip")
    push!(all_data, avg_energies)
    plot!(p, beta_values, avg_energies, label = "Single flip")

    
    xlabel!(p,"Beta", xlabelcolor = :white)
    ylabel!(p,"Energy", ylabelcolor = :white)
    title!(p,"N = $N, copies = $copies, $move, k = $k", titlecolor = :white)
    savefig(joinpath(folder, filename))
    
    # Save data to file
    writedlm(joinpath(folder, datafile), hcat(beta_values, all_data...), ',')
end

N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.5, 100)
k_values = [1,2,3,4,5,10,11,25,26,50,51,75,99]
copies = 10
folder = ""
monte_carlo_timesteps = 10

move = "k chain flip"
n_multiplier = [0.1,0.5, 1.0, 2.0, 4.0, 8.0]
decorrelation_copies = 10
println("started")
for k in k_values
    file_name =  "tau k$k chain flip.png"
    data_file = "tau k$k chain flip.csv"


    time = now()
    println("k=$k")
    figure_E_vary_tau(lattice, beta_values, k, copies, file_name, data_file, folder, N, n_multiplier, move, decorrelation_copies, monte_carlo_timesteps)
    println(time - now())
end
println("Finished all")
