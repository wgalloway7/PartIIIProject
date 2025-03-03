using Plots
using DelimitedFiles
include("lattice.jl")
include("montecarlo.jl")
N = 20
res = 100
#beta_values= 1 ./ generate_T_intervals(10.0,0.2,20)
beta_values = 1 ./ collect(range(0.1, stop=2.0, length=res))

data = readdlm("autocorrelation.csv", ',')
data = replace(data, "" => NaN)
tau_values = []
println(size(data))
for i in 1:size(data)[1]
    row = filter(!isnan, data[i, :])
    print(i)
    p = plot()
    plot!(p, row)
    title!(p, "beta = $(beta_values[i])")
    savefig(p, "autocorrelation$i.png")
    tau = findfirst(x -> x < exp(-1), row)
    push!(tau_values, something(tau, 0))
end
println(tau_values)
println(beta_values)
println(size(tau_values))
println(size(beta_values))
p = plot()
plot!(p, 1 ./ beta_values, tau_values)
vline!(p, [1 / log(1 + sqrt(2))], label = "Tc")
savefig(p, "tau_values.png")