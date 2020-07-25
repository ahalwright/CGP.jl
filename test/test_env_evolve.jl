# Simple test of function run_mut_evolution()
# include("../src/CGP.jl") # Assumes that CGP.jl is loaded
# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
date = "7_25"
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/data/$date/test_env_evolveA.csv"
else
  csvfile = "../data/$date/test_env_evolveA.csv"
end
iterations = 2
numinputs = 3
numoutputs = 6 
nodearity = 2
numinteriors = 15
ngoals = 1
levelsback=15
#hamming_rng = true:true
maxsteps = 100000:100000
gl_repetitions = 3
num_flip_bits = 1:2
perm_heuristic=true
perturb_goal_range = false:true
avgfitness = true
#run_env_evolution( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback, csvfile )
df = run_env_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, 
    gl_repetitions, num_flip_bits, avgfitness, perm_heuristic, perturb_goal_range, csvfile )
