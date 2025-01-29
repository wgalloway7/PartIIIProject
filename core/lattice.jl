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

#calculate energy
#generate matrix of sum of nearest neighbours of each site, multiply by site value to get energy contribution, sum
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


function single_flip!(lattice::Lattice;reverse::Bool=false, candidate_reversing_information=nothing,k::Int64 = 1)
    if !reverse
        x = rand(1:lattice.N)
        y = rand(1:lattice.N)
        lattice.grid[x,y] *= -1
        return (x,y)
    else
        x,y = candidate_reversing_information
        lattice.grid[x,y] *= -1
        #thankfully to get inverse we just flip again
    end
end

function k_chain_flip!(lattice::Lattice;reverse::Bool=false, candidate_reversing_information=nothing,k::Int64=1)
    if !reverse
        x = rand(1:lattice.N)
        y = rand(1:lattice.N)
        #wlog chain is on x axis
        for i in 0:k-1
            lattice.grid[mod1(x+i,N),y] *= -1
        end
        return (x,y)
    else
        x,y = candidate_reversing_information
        for i in 0:k-1
            lattice.grid[mod1(x+i,N),y] *= -1
        end
    end
end

function k_line_flip!(lattice::Lattice;reverse::Bool=false, candidate_reversing_information=nothing,k::Int64=1)
    if !reverse
        y = rand(1:lattice.N)
        xvals=  rand(1:lattice.N,k)
        for x in xvals
            lattice.grid[x,y] *= -1
        end
        return (xvals,y)
    else
        xvals,y = candidate_reversing_information
        for x in xvals
            lattice.grid[x,y] *= -1
        end
    end
end

function unconstrained_line_flip!(lattice::Lattice; reverse::Bool=false, candidate_reversing_information=nothing,k::Int64 =1)
    if !reverse
        xvals = rand(1:lattice.N,k)
        yvals = rand(1:lattice.N,k)
        for i in 1:k
            lattice.grid[xvals[i],yvals[i]] *= -1
        end
        return (xvals,yaxis)
    else
        xvals, yvals = candidating_reversing_information
        for i in 1:k
            lattice.grid[xvals[i],yvals[i]] *= -1
        end
    end
end

function k_square_flip!(lattice::Lattice; reverse::Bool=false, candidate_reversing_information=nothing, k::Int64=1)
    if !reverse
        has_hole = rand(Bool)
        x = rand(1:lattice.N)
        y = rand(1:lattice.N)
        xhole = rand(1:lattice.N)
        yhole = rand(1:lattice.N)
        for i in 0:k-1
            for j in 0:k-1
                if !(x+i == xhole && y+j == yhole && has_hole)
                    lattice.grid[mod1(x+i, lattice.N), mod1(y+j, lattice.N)] *= -1
                end
            end
        end
        return (x, y, xhole, yhole, has_hole)
    else
        x, y, xhole, yhole, has_hole = candidate_reversing_information
        for i in 0:k-1
            for j in 0:k-1
                if !(x+i == xhole && y+j == yhole && has_hole)
                    lattice.grid[mod1(x+i, lattice.N), mod1(y+j, lattice.N)] *= -1
                end
            end
        end
    end
end