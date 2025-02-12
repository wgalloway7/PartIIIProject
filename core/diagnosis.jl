using Random
using Plots
using Statistics
using DelimitedFiles
using Pkg
using LsqFit
include("lattice.jl")
include("montecarlo.jl")
include("main.jl")
using Profile, ProfileView
Profile.view()

