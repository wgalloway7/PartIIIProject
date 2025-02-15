# singleflip.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates
using ForwardDiff
using QuadGK

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
    #C_fluctuations = calculate_C_from_fluctuations(energy_runs, beta_values)
    
    p = plot()
    plot!(p, background_color = "#333333", gridcolor = :white, legend = :topright)
    plot!(p, beta_values[1:end-1], abs.(C_derivative), label = "Derivative", lw = 2)
    #plot!(p, beta_values, abs.(C_fluctuations * N), label = "Fluctuations", lw = 2)
    
    xlabel!(p, "Beta", xlabelcolor = :white)
    ylabel!(p, "C", ylabelcolor = :white)
    title!(p, "Heat capacity from single flip energy runs, N = $N", titlecolor = :white)
    savefig(joinpath(folder, filename))
end
N = 200
beta_values = 1 ./ generate_T_intervals(10.0, 0.25, 100)
#E_runs = readdlm("experiments\\singleflip\\single_flips.csv", ',', Float64)
E_runs = readdlm("experiments\\singleflip\\single_flips_N200.csv", ',', Float64)
E_runs_vector = [collect(row) for row in eachrow(E_runs)]
figure_C(E_runs_vector, beta_values, "C_singleflip.png", "experiments\\singleflip", N)


function entropy_from_heat_capacity(C::Vector{Float64}, beta_values::Vector{Float64})
    beta_values = beta_values[1:end-1]
    T_values = 1.0 ./ beta_values

    reverse!(C)
    reverse!(T_values)

    S = Float64[0.0]
    for i in 2:length(T_values)
        dT = T_values[i] - T_values[i-1]  # Temperature difference
        avg_C_over_T = (C[i]/T_values[i] + C[i-1]/T_values[i-1]) / 2  # Trapezoidal rule
        push!(S, S[end] + avg_C_over_T * dT)  # Integrate (C/T) dT
    end
    return S, T_values
end


# Onsager free energy per spin (with J = 1, k_B = 1)
function onsager_free_energy(b, N)
    # Here K = b (since J = 1)
    K = b
    # Define κ = 2 sinh(2K) / cosh(2K)^2
    kappa = 2 * sinh(2*K) / (cosh(2*K)^2)
    
    # Define the integrand:
    # ln[(1 + sqrt(1 - κ² sin²θ)) / 2]
    function integrand(theta)
        return log(1 + sqrt(1 - kappa^2 * sin(theta)^2))
    end
    
    # Compute the integral from 0 to π and include the 1/(2π) factor
    I, _ = quadgk(integrand, 0, pi/2)
    I /= pi
    
    # Standard form of Onsager's free energy:
    # -b * f = ln(2 cosh(2b)) + I
    # So, f = -1/b * [ln(2 cosh(2b)) + I]
    return (log(sqrt(2) * cosh(2*K)) + I) / (-b)
end

# Onsager energy per spin: E = d/db (b * f)
function onsager_energy(b, N)
    E = ForwardDiff.derivative(b -> b * onsager_free_energy(b, N), b)
    return E
end

# Onsager specific heat per spin: C = b² * dE/db
function onsager_specific_heat(b, N)
    dE_db = ForwardDiff.derivative(b -> onsager_energy(b, N), b)
    C = -b^2 * dE_db
    return C
end

# Onsager entropy per spin: S = b * (E - f)
function onsager_entropy(b, N)
    U = onsager_energy(b, N)
    F = onsager_free_energy(b, N)
    S = b * (U - F)
    return S
end



C = calculate_C_from_derivative(E_runs_vector, beta_values)
S_values,T_vals = entropy_from_heat_capacity(C, beta_values)
energies = mean(hcat(E_runs_vector...), dims=2)[:]


#plot2 = plot()
#plot!(plot2, T_vals, S_values)  # Plot entropy vs. temperature
#xlabel!(plot2, "T", xlabelcolor = :white)
#ylabel!(plot2, "S(T)", ylabelcolor = :white)
#title!(plot2, "Entropy from heat capacity, N = $N", titlecolor = :white)
#savefig(plot2,"experiments\\singleflip\\entropy.png")


plot3 = plot()
#plot!(plot3, 1 ./beta_values, energies, label = "Single flip energy", lw = 2)
plot!(plot3, beta_values, onsager_free_energy.(beta_values,N), label = "Onsager free energy", lw = 2)
xlabel!(plot3, "Beta", xlabelcolor = :white)
ylabel!(plot3, "F(T)", ylabelcolor = :white)
title!(plot3, "onsager free energy, N = $N", titlecolor = :white)
savefig(plot3,"experiments\\singleflip\\onsager free energy.png")


plot4 = plot()
plot!(plot4, beta_values, onsager_energy.(beta_values,N), label = "Onsager energy", lw = 2)
xlabel!(plot4, "Beta", xlabelcolor = :white)
ylabel!(plot4, "E(T)", ylabelcolor = :white)
title!(plot4, "onsager energy, N = $N", titlecolor = :white)
savefig(plot4,"experiments\\singleflip\\onsager energy.png")

plot5 = plot()
plot!(plot5, beta_values, onsager_specific_heat.(beta_values,N), label = "Onsager specific heat capacity", lw = 2)
plot!(plot5, beta_values[1:end-1], C, label = "Specific heat capacity from single flip", lw = 2)
xlabel!(plot5, "Beta", xlabelcolor = :white)
ylabel!(plot5, "C(T)", ylabelcolor = :white)
title!(plot5, "onsager heat capacity, N = $N", titlecolor = :white)
savefig(plot5,"experiments\\singleflip\\onsager heat capacity.png")

plot6 = plot()
plot!(plot6, beta_values, onsager_entropy.(beta_values,N), label = "Onsager entropy", lw = 2)
plot!(plot6, 1 ./T_vals, S_values, label = "Entropy from heat capacity", lw = 2)
xlabel!(plot6, "Beta", xlabelcolor = :white)
ylabel!(plot6, "S(T)", ylabelcolor = :white)
title!(plot6, "onsager entropy, N = $N", titlecolor = :white)
savefig(plot6,"experiments\\singleflip\\onsager entropy.png")

plot7 = plot()
plot!(plot7, 1 ./beta_values, energies, label = "Single flip energy", lw = 2)
plot!(plot7, 1 ./beta_values, onsager_energy.(beta_values,N), label = "Onsager energy", lw = 2)
xlabel!(plot7, "T", xlabelcolor = :white)
ylabel!(plot7, "E(T)", ylabelcolor = :white)
title!(plot7, "onsager energy, N = $N", titlecolor = :white)
savefig(plot7,"experiments\\singleflip\\single_flip_energy.png")