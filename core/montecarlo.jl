#montecarlo.jl
include("lattice.jl")

function monte_carlo_timestep!(lattice::Lattice, move::String, beta::Float64; verbose::Bool=false, energy_cap::Float64=0.0, k::Int64=1)
    flips = generate_moves(lattice, move, k)

    energy_difference = energy_change(lattice,flips)

    if energy_difference <= 0 || rand() < exp(-beta * energy_difference)
        do_flips(lattice, flips)
        return 1
    else
        return 0
    end

end

function run_metropolis_algorithm(lattice::Lattice, beta::Float64, k::Int64, move::String = "single flip"; maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1), verbose::Bool=false)
    #need to get rid of candidate_generating_function
    current_iteration = 0
    accepted_candidates = 0

    initial_lattice = Lattice(lattice.N)
    initial_lattice.grid = deepcopy(lattice.grid)
    current_configuration_correlation_function_value = configuration_correlation_function(lattice, initial_lattice)

    while (current_iteration <= maximum_iterations) && (current_configuration_correlation_function_value > configuration_correlation_convergence_criteria)
 
        accepted_candidates_increase = monte_carlo_timestep!(lattice, move, beta; verbose=verbose, k=k)

        current_iteration += 1
        accepted_candidates += accepted_candidates_increase
        current_configuration_correlation_function_value = configuration_correlation_function(lattice, initial_lattice)
    end

    converged = (current_iteration < maximum_iterations)
    if verbose
        if converged
            printstyled("Configuration correlation function converged!", color = :green)
        else
            printstyled("Maximum Iterations Reached", color = :red)
        end
    end

    return (converged, current_configuration_correlation_function_value, current_iteration, accepted_candidates)
end