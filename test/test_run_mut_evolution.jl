# Simple test of function run_mut_evolution()
# include("../src/CGP.jl") # Assumes that CGP.jl is loaded
# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_run_mut_evolve.csv"
else
  csvfile = "../test/data/test_run_mut_evolve.csv"
end
iterations = 4
numinputs = 2:2
numoutputs = 2:2
nodearity = 2
numinteriors = 6:6
ngoals = 8:8
levelsback=6:6
hamming_rng = true:true
maxsteps = 100000:100000
gl_repetitions=2
#run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback, hamming_rng, csvfile )
run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, hamming_rng, csvfile, gl_repetitions=gl_repetitions )