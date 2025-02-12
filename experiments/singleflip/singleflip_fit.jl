# singleflip.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function calculate_C_from_derivative(energy_runs::Vector{Vector{Float64}}, beta_values::Vector{Float64})

    energies = mean(hcat(energy_runs...), dims=2)[:]
    dE_dBeta = Float64.(diff(energies) ./ diff(beta_values))
    beta_mid = Float64.((beta_values[1:end-1] .+ beta_values[2:end]) ./ 2)
    C = - dE_dBeta .* beta_mid.^2
    return C
end


function calculate_C_from_fluctuations(energy_runs::Vector{Vector{Float64}}, beta_values::Vector{Float64})
    # Initialize a vector to store the heat capacity for each beta value
    C_values = Float64[]

    # Iterate over each beta value
    for i in 1:length(beta_values)
        beta = beta_values[i]
        
        # Collect the energy values for the current beta value (from all runs)
        energies_at_beta = [energy_runs[j][i] for j in 1:length(energy_runs)]
        
        # Calculate the mean and mean squared energy
        E_avg = mean(energies_at_beta)
        E_squared_avg = mean(energies_at_beta .^ 2)
        
        # Calculate the fluctuation of energy (heat capacity)
        C = (E_squared_avg - E_avg^2) * beta^2
        
        # Store the calculated heat capacity
        push!(C_values, C)
    end

    # Return the vector of heat capacity values for each beta
    return C_values
end

function figure_C(energy_runs::Vector{Vector{Float64}}, beta_values::Vector{Float64}, filename::String, folder::String, N::Int64)
    C_derivative = calculate_C_from_derivative(energy_runs, beta_values)
    C_fluctuations = calculate_C_from_fluctuations(energy_runs, beta_values)
    
    p = plot()
    plot!(p, background_color = "#333333", gridcolor = :white, legend = :topright)
    plot!(p, beta_values[1:end-1], abs.(C_derivative), label = "Derivative", lw = 2)
    plot!(p, beta_values, abs.(C_fluctuations * N), label = "Fluctuations", lw = 2)
    
    xlabel!(p, "Beta", xlabelcolor = :white)
    ylabel!(p, "C", ylabelcolor = :white)
    title!(p, "Heat capacity from single flip energy runs, N = $N, copies = $copies", titlecolor = :white)
    savefig(joinpath(folder, filename))
end

beta_values = 1 ./ generate_T_intervals(10.0, 0.25, 100)
E_runs = readdlm("experiments\\singleflip\\single_flips.csv", ',', Float64)
E_runs_vector = [collect(row) for row in eachrow(E_runs)]
println(typeof(E_runs_vector))
figure_C(E_runs_vector, beta_values, "C_singleflip.png", "experiments\\singleflip", N)


