# Simple test of function run_mut_evolution()
# include("../src/CGP.jl") # Assumes that CGP.jl is loaded
# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_run_mut_evolve.csv"
else
  csvfile = "../test/data/test_run_mut_evolve.csv"
end
iterations = 100
numinputs = 4
numoutputs = 1
nodearity = 2
numinteriors = 16
ngoals = 1
levelsback=numinteriors-4
hamming_rng = true
fault_tol_rng = false
maxsteps = 100000
gl_repetitions=1
#fit_limit_list = [collect(numoutputs)[1]-2.0, collect(numoutputs)[1]-1.0, Float64(collect(numoutputs)[1])]
fit_limit_list = [Float64(collect(numoutputs)[1])]
#run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback, hamming_rng, fit_limit_list::Vector{Float64}, csvfile )
df =run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, hamming_rng, fit_limit_list, csvfile, gl_repetitions=gl_repetitions, fault_tol_rng=fault_tol_rng )
