#montecarlo.jl
include("lattice.jl")

function monte_carlo_timestep!(lattice::Lattice, candidate_generating_function!::Function, beta::Float64; verbose::Bool=false, energy_cap::Float64=0.0, k::Int64=1)
    current_energy = energy(lattice)

    if verbose
        println("Current configuration:")
        println(lattice.grid)
    end

    candidate_reversing_information = candidate_generating_function!(lattice; k=k)
    candidate_energy = energy(lattice)

    if verbose
        println("Candidate configuration:")
        println(lattice.grid)
    end

    energy_difference = candidate_energy - current_energy

    if energy_difference <= 0 || rand() < exp(-beta * energy_difference)
        return 1
    else
        candidate_generating_function!(lattice; reverse=true, candidate_reversing_information=candidate_reversing_information, k=k)
        return 0
    end
end

function run_metropolis_algorithm(lattice::Lattice, beta::Float64, k::Int64, move::String = "single flip"; maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1), verbose::Bool=false)
    current_iteration = 0
    accepted_candidates = 0

    initial_lattice = Lattice(lattice.N)
    initial_lattice.grid = deepcopy(lattice.grid)
    current_configuration_correlation_function_value = configuration_correlation_function(lattice, initial_lattice)

    while (current_iteration <= maximum_iterations) && (current_configuration_correlation_function_value > configuration_correlation_convergence_criteria)
        if move == "single flip"
            candidate_generating_function! = single_flip!
        elseif move == "k chain flip"
            candidate_generating_function! = k_chain_flip!
        elseif move == "k line flip"
            candidate_generating_function! = k_line_flip!
        elseif move == "unconstrained k flip"
            candidate_generating_function! = unconstrained_flip!
        elseif move == "k square flip"
            candidate_generating_function! = k_square_flip!
        end

        accepted_candidates_increase = monte_carlo_timestep!(lattice, candidate_generating_function!, beta; verbose=verbose, k=k)

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