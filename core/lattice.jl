#lattice.jl
using Random
using LinearAlgebra
using Distributions 

mutable struct Lattice
    #configuration of Lattice (2D array of binary spins)
    grid::Matrix{Int64}

    N:: Int64

    Lattice(N) = new(solved_configuration(N),N)
end


function solved_configuration(N::Int64)
    return fill(1,N,N)
end

#m average magnetisation [0,1]
function random_configuration(N::Int, m::Float64)
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

function configuration_correlation_function(lattice_configuration::Matrix{Int64}, reference_lattice_configuration::Matrix{Int64})
    if size(lattice_configuration) != size(reference_lattice_configuration)
        throw(ArgumentError("Reference lattice has different size from current lattice"))
    end

    match_count = sum(lattice_configuration .== reference_lattice_configuration)
    total_sites = prod(size(lattice_configuration))
    

    return match_count/ total_sites
end

function energy(lattice::Lattice)
    N = lattice.N
    grid = lattice.grid
    total_energy = 0
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
    return total_energy
end


function generate_moves(lattice::Lattice, move::String, k::Int64 = 1)
    N  = lattice.N
    if move == "single flip"
        x = rand(1:lattice.N)
        y = rand(1:lattice.N)
        return ([x],[y])

    elseif move == "k chain flip"
        x_0 = rand(1:lattice.N)
        y_0 = rand(1:lattice.N)
        println(([mod1(x_0 + i, N) for i in 0:k-1], [y_0 for i in 0:k-1]))
        return ([mod1(x_0 + i, N) for i in 0:k-1], [y_0 for i in 0:k-1])


    elseif move == "k line flip"
        y = rand(1:lattice.N)
        return (rand(1:lattice.N,k), [y for i in 1:k],)
      
    elseif move == "unconstrained k flip"
        return (rand(1:lattice.N,k), rand(1:lattice.N,k))
        
    elseif move == "k square flip"
        has_hole = rand(Bool)
        x0 = rand(1:lattice.N)
        y0 = rand(1:lattice.N)
        
        indices = [(mod1(x0 + i, lattice.N), mod1(y0 + j, lattice.N)) for i in 0:k-1, j in 0:k-1]

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
    ref_grid = deepcopy(lattice.grid)
    xvals, yvals = flips
    if length(xvals) != length(yvals)
        throw(ArgumentError("xvals and yvals must have the same length"))
    end

    ref_grid[xvals, yvals] .= 0
    total_energy_change = 0
    for (x, y) in zip(eachindex(xvals), eachindex(yvals))
        xval = xvals[x]
        yval = yvals[y]
        
        neighbours = sum([
            ref_grid[mod1(xval+1, N), yval],  # Right neighbor
            ref_grid[mod1(xval-1, N), yval],  # Left neighbor
            ref_grid[xval, mod1(yval+1, N)],  # Top neighbor
            ref_grid[xval, mod1(yval-1, N)]   # Bottom neighbor
        ])
        total_energy_change += -1 * lattice.grid[xval,yval] * neighbours
    end
    return total_energy_change
end

function do_flips(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    xvals, yvals = flips
    for (x,y) in zip(xvals, yvals)
        lattice.grid[x,y] *= -1
    end
end

function magnetisation(lattice::Lattice)
    return sum(lattice.grid) / (lattice.N * lattice.N)
end