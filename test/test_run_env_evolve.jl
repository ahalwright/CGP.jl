# Simple test of function run_mut_evolution()
# include("../src/CGP.jl") # Assumes that CGP.jl is loaded
# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
#include("../src/Env_evolution.jl")
#include("../src/Assignment.jl")
cwd = pwd()
date = "7_25"
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/data/$date/test_env_evolveA.csv"
else
  csvfile = "../data/$date/test_env_evolveA.csv"
end
iterations = 1
numinputs = 3
numoutputs = 6 
nodearity = 2
numinteriors = 15
ngoals = 1
levelsback=15
#hamming_rng = true:true
maxsteps = 100000:100000
gl_repetitions = 3
num_flip_bits = 2:2
perm_heuristic=true
perturb_goal_range = false:true
avgfitness = true
fit_limit_list = [numoutputs-2.0, numoutputs-1.0, Float64(numoutputs)]
#run_env_evolution( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback, csvfile )
df = run_env_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, 
    gl_repetitions, num_flip_bits, avgfitness, perm_heuristic, perturb_goal_range, fit_limit_list, csvfile )
