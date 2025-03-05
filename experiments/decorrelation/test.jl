using Plots
using DelimitedFiles
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
N = 100
res = 100
copies = 5
#beta_values= 1 ./ generate_T_intervals(10.0,0.2,20)
beta_values = 1 ./ collect(range(1.5, stop=3.0, length=res))
data = []
tau_values = zeros(res)
for i in 1:copies
    run = readdlm("autocorrelation_$i.csv", ',')
    for j in 1:res
        row = run[j, :]
        tau = length(row)
        for (m,a) in enumerate(row)
            if a < exp(-1)
                tau = m
                break
            end
        end
        tau_values[j] += tau
    end
end
tau_values ./= copies
println(tau_values)



p = plot()
plot!(p, 1 ./ beta_values, tau_values ./ 10)
vline!(p, [2 / log(1 + sqrt(2))], label = "Tc")
xlabel!("Temperature")
ylabel!("Tau (MC steps)")
title!(p,"Tau vs Temperature for N = $N, copies = $copies")
savefig(p, "tau_values_5runs.png")

q = plot()
plot!(q, beta_values, tau_values ./ 10)
vline!(q, [log(1 + sqrt(2)) / 2], label = "Tc")
xlabel!("Beta")
ylabel!("Tau (MC steps)")
title!(q,"Tau vs Beta for N = $N, copies = $copies")

savefig(q, "tau_values_beta_5runs.png")
