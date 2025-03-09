#lattice.jl
using Random
using LinearAlgebra
using Distributions
using Combinatorics
using DelimitedFiles

mutable struct Lattice
    # Configuration of the lattice (2D array of binary spins)
    grid::Matrix{Int64}
    N::Int64
    tau_values::Vector{Int64}  # Vector to store the tau values read from a CSV

    # Constructor to initialize the lattice with solved configuration and tau values
    function Lattice(N::Int64, csv_file::String = "../../core/tau_values_flat.csv")
        grid = solved_configuration(N)
        tau_values = read_tau_values(csv_file)
        new(grid, N, tau_values)
    end
end

@inline function solved_configuration(N::Int64)
    # minimal energy state
    # all spins aligned
    # acessilbe assuming ergodic
    return fill(1,N,N)
end

function read_tau_values(csv_file::String)
    # Read the CSV file
    data =readdlm(csv_file, ',')
    return Int64.(data[2, :])  
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


function configuration_correlation_function(lattice::Lattice, reference_lattice::Lattice, C_0::Int64)
    # OlD, wrong implementation
    # calculate correlation function
    # C(t) = <s(t) s(0)> - <s(t)> <s(0)>
    # C(t) = 1/N^2 [sum(s_i(t) s_i(0)) - sum(s_i(t)) * sum(s_i(0))]
    lattice_configuration = lattice.grid
    reference_lattice_configuration = reference_lattice.grid

    C_t_C_0 = sum(lattice_configuration .* reference_lattice_configuration)
    C_t = sum(lattice_configuration)
    

    
    return (C_t_C_0 / lattice.N) - (C_t * C_0 / lattice.N^2)
end

function energy(lattice::Lattice)
    # for each lattice sites
    # calculate sum of nearest neighbours
    # by shifting lattice and summing
    N = lattice.N
    grid = lattice.grid
    total_energy = 0
    # vectorised energy calculation
    right = circshift(grid, (1, 0))
    left = circshift(grid, (-1, 0))
    up = circshift(grid, (0, 1))
    down = circshift(grid, (0, -1))
    total_energy = sum(grid .* (right + left + up + down))
    return total_energy / (2 * N^2)
end


function generate_moves(lattice::Lattice, move::String, k::Int64 = 1)
    #TO DO:
    #STOP duplicates being generated
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
        return (sample(1:N, k, replace=false), fill(y, k))
      
    elseif move == "unconstrained k flip"
        # randomly select k sites
        x_sites = sample(1:N, k, replace=false)
        y_sites = sample(1:N, k, replace=false)
        return (x_sites, y_sites)
        
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

@inline function energy_change(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    N = lattice.N
    
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
    flip_set = Set(zip(xvals,yvals))
    total_energy_change = 0
    @inbounds for (x, y) in zip(xvals, yvals)
        local sum_neighbours = 0
        for (dx,dy) in ((1,0),(-1,0),(0,1),(0,-1))
            xp = mod1(x+dx, N)
            yp = mod1(y+dy, N)
            if (xp,yp) ∉ flip_set
                sum_neighbours += lattice.grid[xp,yp]
            end
        end
        total_energy_change += -lattice.grid[x, y] * 2 * sum_neighbours
    end
    # / N completely arbitrary, just shifts temperature scale
    return total_energy_change
end

function do_flips(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    # given list of candidate sites
    # flip the spins of these sites
    # same site might be flipped multiple times?
    xvals, yvals = flips
    @inbounds for (x,y) in zip(xvals, yvals)
        lattice.grid[x,y] *= -1
    end
end

function magnetisation(lattice::Lattice)
    return sum(lattice.grid) / (lattice.N * lattice.N)
end


function explore_moves(lattice::Lattice, k::Int64, move::String)
    moves = []

    if move == "single flip"
        moves = [([i], [j]) for i in 1:lattice.N, j in 1:lattice.N]   

    elseif move == "k chain flip"
        for y in 1:lattice.N
            for x in 1:(lattice.N - k + 1)
                x_coords = [(x + i - 1) % lattice.N + 1 for i in 1:k]
                y_coords = fill(y, k)
                push!(moves, (x_coords, y_coords))
            end
        end

    elseif move == "k line flip"
        for y in 1:lattice.N
            for comb in collect(combinations(1:lattice.N, k))  # Avoid lazy evaluation issues
                push!(moves, (collect(comb), fill(y, k)))
            end
        end
    
    elseif move == "unconstrained k flip"
        indices = [(i, j) for i in 1:lattice.N for j in 1:lattice.N]  # Flatten grid
        for comb in collect(combinations(indices, k))
            x_coords, y_coords = unzip(comb)  # Efficient tuple splitting
            push!(moves, (collect(x_coords), collect(y_coords)))
        end



    else
        throw(ArgumentError("Invalid move type"))
        return 0
    end


    
    # if energy change of proposed move is negative, add -1 to up_down
    # if energy change is positive, add 1 to down_up
    saddles = count(move -> energy_change(lattice, move) < 0, moves)
    return saddles
    
end










