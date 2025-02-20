# VTM.jl
# run to generate fits for decorrelation curves for various beta

using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function fit_VTM(correlation_histories::Matrix{Float64}, beta_values::Vector{Float64}, N::Int64, move::String, filename::String, datafile::String, k::Int64 = 1)
    
    function VTM(t, p)
        tau, alpha = p
        return exp.(- (t ./ tau) .^ alpha)
    end

    params = []
    p = plot()
    colors = cgrad(:RdBu, length(beta_values))
    for (i, beta) in enumerate(beta_values)
        println(i)
        y_data = correlation_histories[i, :]  # Use row i for each beta
        x_data = Float64.(1:length(y_data))

        initial_params = [Float64(N^2), 1.0]

        lower_bounds = [1.0, 1e-8]
        upper_bounds = [Inf, Inf]
        fit = curve_fit(VTM, x_data, y_data, initial_params, lower = lower_bounds, upper = upper_bounds)
        fitted_params = fit.param
        push!(params, vcat(fitted_params, beta))
        tau = round(fitted_params[1])
        alpha = round(fitted_params[2], digits = 2)
        beta_round = round(beta, digits = 2)

        plot!(p, x_data, y_data, label="$beta_round", color = colors[i])
        plot!(p, x_data, VTM(x_data, fitted_params), label = "fit, tau = $tau, alpha = $alpha", linestyle=:dash, color = "blue")

    end
    plot!(p, background_color="#333333", gridcolor=:white, legend=:topright)
    xlabel!(p, "Iterations (t)")
    ylabel!(p, "Autocorrelation C(t)")
    title!(p, "VTM fit for $move, k = $k, N = $N")

    savefig(p, filename)
    writedlm(datafile, hcat(params...),',')
end

function parse_csv(filename::String)
    data_matrix = (readdlm(filename, ','))
    return [parse(Float64, strip(replace(string(data_matrix[i, j]), r"[\[\]\"\s]+" => ""))) for i in 1:size(data_matrix, 1), j in 1:size(data_matrix, 2)]
end

N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)


beta_values = 1 ./ generate_T_intervals(10.0, 0.8, 10)
copies = 50
measurement_MC_steps = 10
cooling_MC_steps = 10

k = 98
move = "k chain flip"

filename =  "VTM_fit_kchain98.png"
datafile =  "VTM_fit_kchain98.csv"

data_matrix = parse_csv("VTM_kchain98.csv")

println("started")
time = now()
fit_VTM(data_matrix, beta_values, N, move, filename, datafile)
println(now() - time)




