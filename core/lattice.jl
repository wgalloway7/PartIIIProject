#lattice.jl
using Random
using LinearAlgebra
using Distributions
using Combinatorics

mutable struct Lattice
    #configuration of Lattice (2D array of binary spins)
    grid::Matrix{Int64}

    N:: Int64

    Lattice(N) = new(solved_configuration(N),N)
end


function solved_configuration(N::Int64)
    # minimal energy state
    # all spins aligned
    # acessilbe assuming ergodic
    return fill(1,N,N)
end

function random_configuration(N::Int, m::Float64)
    # random configuration of spins
    # weighted probability of spin up or down
    if m < -1 || m > 1
        throw(ArgumentError("m out of range"))
    end

    choices = [-1, 1]
    p = 0.5 * (m + 1.0)
    probabilities = [1.0 - p, p]
    spin_distribution = Categorical(probabilities)
    return [rand(spin_distribution) == 2 ? 1 : -1 for _ in 1:N, _ in 1:N]
end


function configuration_correlation_function(lattice::Lattice, reference_lattice::Lattice)
    return configuration_correlation_function(lattice.grid, reference_lattice.grid)
end
#TO FIX: why is there two of these?

function configuration_correlation_function(lattice_configuration::Matrix{Int64}, reference_lattice_configuration::Matrix{Int64})
    # calculates correlation function given lattice and reference lattice


    if size(lattice_configuration) != size(reference_lattice_configuration)
        throw(ArgumentError("Reference lattice has different size from current lattice"))
    end

    match_count = sum(lattice_configuration .== reference_lattice_configuration)
    total_sites = prod(size(lattice_configuration))
    
    #normalising
    return match_count/ total_sites
end

function energy(lattice::Lattice)
    # for each lattice sites
    # calculate sum of nearest neighbours
    # by shifting lattice and summing
    N = lattice.N
    grid = lattice.grid
    total_energy = 0
    # currently O(N^2)
    # could be vectorised by operating shifts directly on lattice
    # larger space complexity but smaller time complexity
    for i in 1:N
        for j in 1:N
            total_energy += grid[i, j] * (
                grid[mod1(i+1, N), j] +  # Right neighbor
                grid[mod1(i-1, N), j] +  # Left neighbor
                grid[i, mod1(j+1, N)] +  # Top neighbor
                grid[i, mod1(j-1, N)]    # Bottom neighbor
            )
        end
    end
    return total_energy / (2 * N^2)
end


function generate_moves(lattice::Lattice, move::String, k::Int64 = 1)
    # generate candidate site for flip moves
    # based on move type
    # returns tuple of x and y coordinates
    N  = lattice.N
    if move == "single flip"
        # randomly select a site
        x = rand(1:lattice.N)
        y = rand(1:lattice.N)
        return ([x],[y])

    elseif move == "k chain flip"
        # randomly select an initial site
        x_0 = rand(1:lattice.N)
        y_0 = rand(1:lattice.N)
        # pick subsequent sites on x axis
        return ([mod1(x_0 + i, N) for i in 0:k-1], [y_0 for i in 0:k-1])


    elseif move == "k line flip"
        # randomly select a line
        y = rand(1:lattice.N)
        # within line randomly select k sites
        return (rand(1:lattice.N,k), [y for i in 1:k],)
      
    elseif move == "unconstrained k flip"
        # randomly select k sites
        return (rand(1:lattice.N,k), rand(1:lattice.N,k))
        
    elseif move == "k square flip"
        # randomly select an initial site
        x0 = rand(1:lattice.N)
        y0 = rand(1:lattice.N)

        # select whether flip has a hole (k^2 - 1 flips) or not (k^2 flips)
        has_hole = rand(Bool)
        
        #generate indices for k x k square
        indices = [(mod1(x0 + i, lattice.N), mod1(y0 + j, lattice.N)) for i in 0:k-1, j in 0:k-1]

        #remove hole if needed
        if has_hole && !isempty(indices)
            hole_index = rand(1:length(indices))
            deleteat!(indices, hole_index)
        end
        xvals = [coord[1] for coord in indices]
        yvals = [coord[2] for coord in indices]
        return (xvals, yvals)
    else
        throw(ArgumentError("Invalid move type"))
    end
end

function energy_change(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    N = lattice.N
    # copy of lattice to mask with proposed flips
    ref_grid = deepcopy(lattice.grid)
    # unpack proposed flips
    xvals, yvals = flips

    if length(xvals) != length(yvals)
        throw(ArgumentError("xvals and yvals must have the same length"))
    end

    # mask proposed flips
    # takes into account the fact that flipping two adjacent spins
    # doesn't change their interaction energy
    # ie +1,+1 has same energy as -1,-1
    # so by masking we acknowledge that the energy change is 0
    ref_grid[xvals, yvals] .= 0
    total_energy_change = 0
    for (x, y) in zip(eachindex(xvals), eachindex(yvals))
        xval = xvals[x]
        yval = yvals[y]
        # calculating sum of neighbours for proposed sites
        # note periodic boundary conditions

        neighbours = sum([
            ref_grid[mod1(xval+1, N), yval],  # Right neighbor
            ref_grid[mod1(xval-1, N), yval],  # Left neighbor
            ref_grid[xval, mod1(yval+1, N)],  # Top neighbor
            ref_grid[xval, mod1(yval-1, N)]   # Bottom neighbor
        ])
        total_energy_change += -1 * lattice.grid[xval,yval] * neighbours
    end
    # / N completely arbitrary, just shifts temperature scale
    return total_energy_change
end

function do_flips(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    # given list of candidate sites
    # flip the spins of these sites
    # same site might be flipped multiple times?
    xvals, yvals = flips
    for (x,y) in zip(xvals, yvals)
        lattice.grid[x,y] *= -1
    end
end

function magnetisation(lattice::Lattice)
    return sum(lattice.grid) / (lattice.N * lattice.N)
end


function explore_moves(lattice::Lattice, k::Int64, move::String)
    moves = []

    if move == "single flip"
        for i in 1:lattice.N
            for j in 1:lattice.N
                push!(moves, ( [i], [j] ))
            end
        end    

    elseif move == "k chain flip"
        for y in 1:lattice.N
            for x in 1:(lattice.N - k + 1)
                x_coords::Vector{Int64} = []
                y_coords::Vector{Int64} = []
                for i in 1:k
                    push!(x_coords, mod1(x + i - 1, lattice.N))
                    push!(y_coords, y)
                end
                push!(moves, (x_coords, y_coords))
            end
        end

    elseif move == "k line flip"
        for y in 1:lattice.N
            for comb in combinations(1:lattice.N, k)
                x_coords::Vector{Int64} = []
                y_coords::Vector{Int64} = []
                for x in comb
                    push!(x_coords, x)
                    push!(y_coords, y)
                end
                push!(moves, (x_coords, y_coords))
            end
        end
    
    elseif move == "unconstrained k flip"
        indices = [(i, j) for i in 1:lattice.N, j in 1:lattice.N]
        for comb in combinations(indices, k)
            x_coords::Vector{Int64} = []
            y_coords::Vector{Int64} = []
            for (i, j) in comb
                push!(x_coords, i)
                push!(y_coords, j)
            end
            push!(moves, (x_coords, y_coords))
        end



    else
        throw(ArgumentError("Invalid move type"))
    end


    
    # if energy change of proposed move is negative, add -1 to up_down
    # if energy change is positive, add 1 to down_up
    saddles = 0
    for move in moves
        move_energy_change = energy_change(lattice, move)
        if move_energy_change < 0
            saddles += 1
        end
    end
    return saddles
    
end










