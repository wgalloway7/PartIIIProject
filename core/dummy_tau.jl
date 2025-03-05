using DelimitedFiles
res = 100
tau_dummy = 200
beta_dummy = 0.0
tau_values = [tau_dummy for _ in 1:res]
beta_values = [beta_dummy for _ in 1:res]
writedlm("tau_values_flat_low.csv", (beta_values, tau_values), ',')