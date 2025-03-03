using Plots
using DelimitedFiles
include("lattice.jl")
include("montecarlo.jl")
N = 20
res = 100
#beta_values= 1 ./ generate_T_intervals(10.0,0.2,20)
beta_values = 1 ./ collect(range(0.2, stop=5.0, length=res))

data = readdlm("autocorrelation_short.csv", ',')
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
plot!(p, 1 ./ beta_values, tau_values ./ 10)
vline!(p, [2 / log(1 + sqrt(2))], label = "Tc")
xlabel!("Temperature")
ylabel!("Tau (MC steps)")
savefig(p, "tau_values.png")
q = plot()
plot!(q, beta_values, tau_values)
vline!(q, [log(1 + sqrt(2)) / 2], label = "Tc")
xlabel!("Beta")
ylabel!("Tau (MC steps)")
savefig(q, "tau_values_beta.png")