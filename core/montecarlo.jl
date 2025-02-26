#montecarlo.jl
include("lattice.jl")

@inline function monte_carlo_timestep!(lattice::Lattice, move::String, beta::Float64; verbose::Bool=false, k::Int64=1)
    # generates potential move's flip sites
    flips = generate_moves(lattice, move, k)

    # calculates energy difference of proposed move
    energy_difference = energy_change(lattice,flips)

    # MH logic on whether to accept the proposed move
    if energy_difference <= 0 || rand() < exp(-beta * energy_difference)
        # if move is accepted, flip the propose sites
        do_flips(lattice, flips)
        return 1
    else
        return 0
    end

end

function run_metropolis_algorithm(lattice::Lattice, beta::Float64, k::Int64, move::String = "single flip"; maximum_iterations::Int64=10000, configuration_correlation_convergence_criteria::Float64=exp(-1.0), verbose::Bool=false, use_correlation::Bool=true, correlation_measurements_per_MC_timstep::Int64=1)
    # runs metropolis-hastings algorithm
    # for given move funciton, beta and k
    # iteration cutoff either predetermined or until correlation function drops to 1/e

    current_iteration = 0
    accepted_candidates = 0

    iterations_per_correlation_measurement = Int64(floor(lattice.N^2 / correlation_measurements_per_MC_timstep))

    initial_lattice = Lattice(lattice.N)
    initial_lattice.grid = copy(lattice.grid)
    #if we are continuously computing the correlation function we need to compute it here
    #if we are using predefined cutoff iteration numbers from previously calculated decorrelations
    #no need to compute it here
    if use_correlation
        current_configuration_correlation_function_value = configuration_correlation_function(lattice, initial_lattice)

    else
        current_configuration_correlation_function_value = 1.0
    end

    #iterate algorithm until cutoff
    while (current_iteration <= maximum_iterations) && (current_configuration_correlation_function_value > configuration_correlation_convergence_criteria)
 
        accepted_candidates_increase = monte_carlo_timestep!(lattice, move, beta; verbose=verbose, k=k)

        current_iteration += 1
        accepted_candidates += accepted_candidates_increase
        #calculate new correlation function after iteration
        # calculate correlation function every lattice.N^2 iterations
        # ie 1 Monte-Carlo iteration
        #if use_correlation
        if use_correlation && current_iteration % iterations_per_correlation_measurement == 0
            current_configuration_correlation_function_value = configuration_correlation_function(lattice, initial_lattice)
        end

        
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

function prepare_lattice!(lattice::Lattice,k::Int64 = 1; maximum_iterations::Int64=10000)
    #prepares lattice in minimal energy state (all aligned)
    #this is an accessible state if system is ergodic
    #then runs metropolis algorithm for beta = 0
    #equivalent to infinite temperature, all moves are accepted
    lattice.grid = solved_configuration(lattice.N)
    run_metropolis_algorithm(lattice, 0.0, k, maximum_iterations = maximum_iterations)
end

function generate_decorrelation_n(lattice::Lattice, beta_values::Vector{Float64}; k::Int64 = 1, move::String = "single flip", maximum_iterations::Int64 = 10000, copies::Int64 = 1)
    decorrelation_n_matrix = zeros(Int64, length(beta_values), copies)

    for copy_idx in 1:copies
        # prepare lattice in a 'hot' state
        # use single flip
        prepare_lattice!(lattice,1 ;maximum_iterations = maximum_iterations)
        
    

        # for each beta value
        # run metropolis algorithm until correlation function drops to 1/e
        # record number of iterations needed
    
        for beta_idx in 1:length(beta_values)
            beta = beta_values[beta_idx]
            (_, _, current_iteration, _) = run_metropolis_algorithm(lattice, beta, k, move, use_correlation = true, maximum_iterations = maximum_iterations)
            
            decorrelation_n_matrix[beta_idx, copy_idx] = current_iteration
        end
    end
    mean_decorrelation_n = mean(decorrelation_n_matrix, dims=2)[:]
    mean_n_int = Int64.(ceil.(mean_decorrelation_n))

    return mean_n_int
    
end

function generate_energies(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, n_correlation::Vector{Int64}, move::String = "single flip")
    #n_correlation = generate_decorrelation_n(lattice, beta_values; k=k, move="single flip")
    energies = Float64[]
    # beta values array of decreasing temperatures
    # as we cool from hot to cold
    # so increasing beta
    
    # for each beta value
    # run metropolis algorithm for n_correlation[i] iterations
    # our decorrelation cutoff determined by "generate_decorrelation_n"
    # calculate the total energy of the lattice and store in array

    #prepare lattice in hot state
    prepare_lattice!(lattice, 1)
    for (i,beta) in enumerate(beta_values)
        #println(beta)
        run_metropolis_algorithm(lattice, beta, k, move, maximum_iterations = n_correlation[i], use_correlation = false)
        E = energy(lattice)
        push!(energies, E)
    end
    return energies
end

function generate_T_intervals(T_hot::Float64, T_cold::Float64, num_intervals::Int64)
    # generates a non-linearly spaced array of temperatures
    # we want a higher density of temperature at T_cold than T_hot
    # as we are interested in the low temperature regime
    if T_hot <= T_cold
        throw(ArgumentError("T_hot must be greater than T_cold"))
    elseif T_cold == 0.0
        throw(ArgumentError("T_cold must be greater than 0"))
    end
    return T_hot .* (T_cold / T_hot) .^ (range(0, 1, length=num_intervals))
end


function generate_correlations(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, move::String = "single flip", max_measurement_iterations::Int64 = 1000, cooling_time::Int64 = 1000)
    # prepare lattice in a 'hot' state
    #TO DO: take average of multiple runs
    prepare_lattice!(lattice, k)
    correlations_beta = []
    # for each beta value
    for beta in beta_values
        # cool to beta
        run_metropolis_algorithm(lattice, beta, k, move, use_correlation = true, maximum_iterations = cooling_time)
        current_iteration = 0
        correlation_history = []
        # making reference lattice to calculate correlation function from
        ref_lattice = Lattice(lattice.N)
        ref_lattice.grid = copy(lattice.grid)
        while current_iteration < max_measurement_iterations
            # run metropolis algorithm for 1000 iterations
            # calculate correlation function and add to history array
            push!(correlation_history, configuration_correlation_function(lattice, ref_lattice))
            monte_carlo_timestep!(lattice, move, beta, k = k)
            current_iteration += 1
        #break after max_measurement_iterations
        end
        push!(correlations_beta, correlation_history)
    end
    return correlations_beta
end

    
function generate_saddles_run(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, move::String = "single flip", cooling_time::Int64 = 10000)
    # prepare lattice in a 'hot state'
    # for each beta
    # cool to beta using single spin flips
    # once cooled, explore moves, calculate energy of system

    prepare_lattice!(lattice)
    # we could also use the decorrelation_n here
    # to determine how many iterations it takes to cool?
    # make sure we're starting hot (low beta)
    beta_values = sort(copy(beta_values))
    saddle_values = Int64[]
    energy_values = Float64[]

    for beta in beta_values
        run_metropolis_algorithm(lattice, beta, k, "single flip", use_correlation = true, maximum_iterations = cooling_time)
        push!(energy_values, energy(lattice))
        saddles = explore_moves(lattice, k, move)
        push!(saddle_values,saddles)
    end
    return (saddle_values, energy_values, beta_values)
end

function generate_saddles(lattice::Lattice, beta_values::Vector{Float64}, k::Int64, move::String,  m::Int64 = 1, cooling_time::Int64 = 10000)
    # generate multiple runs of saddle points
    # and return the average saddle point
    # and energy values
    saddle_values = []
    energy_values = []
    output_beta_vals = []
    critical_energy_values = []
    critical_beta_values = []

   # to use threads?
    for j in 1:m
         #add noise to beta values
        noisy_beta_values = beta_values .+ randn(length(beta_values)) .* 0.1 .* beta_values
        #println(beta_values)
        #println(noisy_beta_values)
        (saddles, energies, betas) = generate_saddles_run(lattice, noisy_beta_values, k, move, cooling_time)
        push!(saddle_values, saddles)
        push!(energy_values, energies)
        push!(output_beta_vals, betas)
        println("run = $j")
        #identifying tempetature and beta at which saddle index vanishes upon cooling.
        critical_energy = -2.0
        critical_beta = 0.0
        for (i,s) in enumerate(saddles)
            if s == 0
                critical_energy = energies[i]
                critical_beta = betas[i]
                break
            end
        end
        push!(critical_energy_values, critical_energy)
        push!(critical_beta_values, critical_beta)
        

    

    
    end
    return (saddle_values, energy_values, output_beta_vals, critical_energy_values, critical_beta_values, k)
end