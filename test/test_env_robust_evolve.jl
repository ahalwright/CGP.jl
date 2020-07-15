# Simple test of function run_mut_evolution()
# include("../src/CGP.jl") # Assumes that CGP.jl is loaded
# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_env_robust_evolve.csv"
else
  csvfile = "../test/data/test_env_robust_evolve.csv"
end
iterations = 40
numinputs = 3
numoutputs = 3 
nodearity = 2
numinteriors = 10
ngoals = 8
levelsback=8
hamming_rng = true:true
maxsteps = 100000:100000
gl_repetitions=2
#run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback, hamming_rng, csvfile )
df = run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, hamming_rng, csvfile, gl_repetitions=gl_repetitions )
