#saddles.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

N = 100
lattice = Lattice(N)
lattice.grid = solved_configuration(N)
beta_values = 1 ./ generate_T_intervals(10.0, 0.25, 100)
folder = ""

k_values = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,45,50,51,55,60,65,70,75,80,85,90,95,99]
#k_values = [98,99]

critical_beta = []
critical_energy = []
move = "k chain flip"

println("started")
for k in k_values
    println("k=$k")
    time = now()
    #datafile = "test.csv"
    datafile = "saddles_k_chain_flip$k.csv"
    results = generate_saddles(lattice, beta_values, k, move,10)
    push!(critical_energy, mean(results[4]))
    push!(critical_beta, mean(results[5]))

    writedlm(joinpath(folder, datafile), results,',')
    println(now() - time)
    p1 = plot()
    xlabel!(p1, "Energy")
    ylabel!(p1, "Average saddle index density")
    title!(p1, "$move, k = $k")
    scatter!(p1, results[2], results[1], markerstrokedwidth = 1, markersize = 1, label = "")
    savefig(p1, joinpath(folder, "saddles_k_chain_flip_energy$k.png"))

    p2 = plot()
    xlabel!(p2, "Beta")
    ylabel!(p2, "Average saddle index density")
    title!(p2, "$move, k = $k")
    scatter!(p2, results[3], results[1], markerstrokedwidth = 1, markersize = 1, label = "")
    savefig(p2, joinpath(folder, "saddles_k_chain_flip_beta$k.png"))
    
end

println("critE = $critical_energy")
println("critb = $critical_beta")
println("k= $k_values")
p3 = plot()
xlabel!(p3, "k")
ylabel!(p3, "Critical Energy")
title!(p3, "Critical Energy vs k")
scatter!(p3, k_values, critical_energy, label = "")
savefig(p3, joinpath(folder, "critical_saddles_k_chain_flip_energy.png"))

p4 = plot()
xlabel!(p4, "k")
ylabel!(p4, "Critical Beta")
title!(p4, "Critical Beta vs k")
scatter!(p4, k_values, critical_beta, label = "")
savefig(p4, joinpath(folder, "critical_saddles_k_chain_flip_beta.png"))

writedlm("critical_saddles.csv",(k_values, critical_energy, critical_beta),',')
println("finished")