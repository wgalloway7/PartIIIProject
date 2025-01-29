#montecarlo.jl
include("lattice.jl")

function monte_carlo_timestep!(lattice::Lattice, candidate_generating_function!::Function, beta::Float64; verbose::Bool=false, energy_cap::Float64=0.0)
    current_energy = energy(lattice)

    if verbose
        println("Current configuration:")
        println(lattice.grid)
    end

    candidate_reversing_information = candidate_generating_function!(lattice)
    candidate_energy = energy(lattice)

    if verbose
        println("Attempted move:")
        println(candidate_reversing_information)
    end


    if verbose
        println("Current Energy: $current_energy")
        println("Candidate Energy: $candidate_energy")
        println("Alpha: $(exp(beta*(current_energy-candidate_energy)))")
    end

    if energy_cap != 0.0 && candidate_energy > energy_cap
        if verbose
            printstyled("Candidate energy above energy cap \n"; color=red)
        end

        candidate_generating_function!(lattice;reverse=true,candidate_reversing_information=candidate_reversing_information)
        
        #return accepted_candidates_increase = 0
        return 0
    end

    if beta == 0.0
        return 1
    end

    
    #calculate acceptance probability, compare to random [0,1]
    alpha = exp(beta * (current_energy - candidate_energy))

    if rand() <= alpha
        if verbose
            printstyled("Switched \n"; color =:green)
            if alpha < 1.0
                printstyled("To higher energy \n"; color = red)
                #not sure what is going on here
            end
        end
        #also return accepted_candidates_increase = 1
        return 1
    else
        #reject candidate and revert it back to the original configuration
        candidate_generating_function!(lattice;reverse=true,candidate_reversing_information=candidate_reversing_information)
        #also return accepted_candidates_increase = 0
        return 0
    end
end



function run_metropolis_algorithm(lattice::Lattice, beta::Float64, maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1); verbose::Bool=false)
    current_iteration = 0
    accepted_candidates = 0

    initial_lattice = Lattice(lattice.N)
    initial_lattice.grid = deepcopy(lattice.grid)
    current_configuration_correlation_function_value = configuration_correlation_function(lattice, initial_lattice)

    while (current_iteration <= maximum_iterations) && (current_configuration_correlation_function_value > configuration_correlation_convergence_criteria)
        move_type = "single flip"
        if move_type == "single flip"
            candidate_generating_function! = single_flip!
        end

        accepted_candidates_increase = monte_carlo_timestep!(lattice, candidate_generating_function!, beta; verbose=verbose)

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
