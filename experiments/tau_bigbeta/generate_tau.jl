using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
using Base.Threads
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


N = 100
res = 200
lattice = Lattice(N)
beta_values= 1 ./ generate_T_intervals(10.0,0.15,res)
#beta_values = 1 ./ collect(range(0.2,stop = 10.0, length=res))
time = now()
filename = "autocorr_b1_5.csv"
equilib_steps = 1000
measurement_steps = 2000
max_lag = 250
copies = 5

println("Starting")
println("Current time is $(now())")



function run_parallel_simulation(lattice::Lattice, beta_values::Vector{Float64}, filename::String, equilib_steps::Int64, measurement_steps::Int64, max_lag::Int64, copies::Int64)
    tau_runs = []
    E_runs = []
    m_runs = []

    # Parallelize the loop over copies
    lock_obj = ReentrantLock()
    @threads for i in 1:copies
        println("Starting copy $(i)")
        time = now()

        # Create a new lattice instance for each thread to ensure no shared state between threads
        thread_lattice = Lattice(lattice.N)  # Ensure each thread uses its own lattice

        # Call the generate_tau function for each copy (parallelized per copy)
        tau_values, E, m = generate_tau(thread_lattice, beta_values, true, filename = "$(i)_$(filename)", max_lag = max_lag, equilib_steps = equilib_steps, measurement_steps = measurement_steps)

        # Store the result for this copy
        lock(lock_obj)
        println("Finished copy $(i) in $(now() - time)")
        try
            push!(tau_runs, tau_values)
            push!(E_runs, E)
            push!(m_runs, m)
        finally
            unlock(lock_obj)
        end
    end

    return tau_runs, E_runs, m_runs
end


tau_runs, E_runs, m_runs = run_parallel_simulation(lattice, beta_values, filename, equilib_steps, measurement_steps, max_lag, copies)
average_tau = mean(hcat(tau_runs...), dims = 2)[:]
average_E = mean(hcat(E_runs...), dims = 2)[:]
average_m = mean(hcat(m_runs...), dims = 2)[:]
println(now() - time)
writedlm("average_tau_b2_5.csv", (beta_values, average_tau, average_E, average_m), ',')

p = plot()
plot!(p, 1 ./ beta_values, average_tau)
plot!(p, 1 ./ beta_values, average_E * -10)
plot!(p, 1 ./ beta_values, average_m * 10)
vline!(p, [2 / log(1 + sqrt(2))], label = "Tc")
vline!(p, [log(1 + sqrt(2)) / 2], label = "1/Tc")

xlabel!("Temperature")
ylabel!("Tau (MC steps)")
savefig(p, "tau_values_b2_5.png")
q = plot()
plot!(q, beta_values, average_tau)
plot!(q, beta_values, average_E * -10)
plot!(q, beta_values, average_m * 10)
vline!(q, [log(1 + sqrt(2)) / 2], label = "Tc")
vline!(q, [2/ log(1 + sqrt(2))], label = "1/Tc")
xlabel!("Beta")
ylabel!("Tau (MC steps)")
savefig(q, "tau_values_beta_b2_5.png")