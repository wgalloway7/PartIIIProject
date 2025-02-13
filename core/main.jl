#main.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Pkg
using LsqFit
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





function figure_correlation_decay(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, filename::String, max_measurement_iterations::Int64, m::Int64, cooling_time::Int64 = 1000, move::String = "single flip", do_plot = true)
    colors = cgrad(:RdBu, length(beta_values))
    
    
    # Initialize storage for correlation histories
    correlation_histories = [ [] for _ in beta_values ]

    # Run multiple instances and collect data
    # m instances
    for _ in 1:m
        correlations_beta = generate_correlations(lattice, beta_values, k, move, max_measurement_iterations, cooling_time)
        for (i, corr_history) in enumerate(correlations_beta)
            push!(correlation_histories[i], corr_history)  # Store runs
        end
    end

    # Compute mean correlation decay for each beta
    averaged_correlations = [ mean(reduce(hcat, runs), dims=2)[:] for runs in correlation_histories ]

    # Plot results
    if do_plot
        p = plot()
        plot!(p, background_color="#333333", gridcolor=:white, legend=:topright)
        for (i, beta) in enumerate(beta_values)
            plot!(p, averaged_correlations[i], label="β = $beta", color=colors[i])
        end
    

        savefig(p, filename)
    end
    return averaged_correlations
end

#N = 50
#lattice = Lattice(N)
#lattice.grid = solved_configuration(N)

#beta_values = 1 ./ generate_T_intervals(10.0, 0.8, 6)
#filename = "hmmm.png"
#m = 10

#figure_correlation_decay(lattice, beta_values, 1, filename, 10000, m, 10000, "single flip")



function fit_VTM(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, filename::String, folder::String, datafile::String, move::String = "single flip", max_measurement_iterations::Int64 = 1000, cooling_time::Int64 = 1000, m::Int64=1, output::Bool = true)
    averaged_correlations = figure_correlation_decay(lattice, beta_values, k, "VTM.png", max_measurement_iterations, m, cooling_time, move, false)
    # Fit the VTM
    # VTM: C(t) = A * exp(-(t/tau)^beta) + C
    # C(t = infinity) = 0
    # C(t = 0) = 1
    # So A = 1, C = 0

    function VTM(t, p)
        tau, alpha = p
        return exp.(-(t ./ tau) .^ alpha)
    end
    p = plot()
    plot!(p, background_color="#333333", gridcolor=:white, legend=:topright)
    colors = cgrad(:RdBu, length(beta_values))
    params = []
    for (i,beta) in enumerate(beta_values)
        y_data = averaged_correlations[i]
        x_data = Float64.(1:length(y_data))

        initial_params = [1000.0, 1.0]

        fit = curve_fit(VTM, x_data, y_data, initial_params)
        fitted_params = fit.param
        push!(params, vcat(fitted_params, beta))
        tau = round(fitted_params[1])
        alpha = round(fitted_params[2], digits=2)
        beta_round = round(beta, digits = 2)

        plot!(p, x_data, y_data, label="beta = $beta_round", linewidth=2, color = colors[i])
        plot!(p, x_data, [VTM(t, fitted_params) for t in x_data], label="fit, tau = $tau, alpha = $alpha", linewidth=2, linestyle=:dash, color = "blue")
        xlabel!(p,"Iterations (t)", xlabelcolor = :white)
        ylabel!(p,"Autocorrelation C(t)", ylabelcolor = :white)
        title!(p,"VTM fit for $move, k = $k", titlecolor = :white)

    end
    savefig(p, joinpath(folder, filename))
    writedlm(joinpath(folder, datafile), hcat(params...), ',')

    return params
end