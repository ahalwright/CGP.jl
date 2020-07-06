# Test the mut_evolve_subgoal() function from src/Evolve.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = Main.CGP.MyInt

# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_inf_alleles.csv"
else
  csvfile = "../test/data/test_inf_alleles.csv"
end
nreps = 3
numinputs = 2
numoutputs = 2
numints = 6
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numints, numlevelsback=numints )
context = construct_context(numinputs)
popsize = 5
max_pop_gens = 50000
indiv_steps = 30
robust_steps = 7
fit_steps = 7
ngoals = 8   # number goals in goallist
tourn_size = 0  
gl = randgoallist(ngoals,p.numinputs,p.numoutputs)
run_robust_evolve( nreps, p, popsize, max_pop_gens, indiv_steps, robust_steps,fit_steps, gl,tourn_size, csvfile )

