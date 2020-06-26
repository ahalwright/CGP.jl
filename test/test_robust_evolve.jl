# Test the mut_evolve_subgoal() function from src/Evolve.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = Main.CGP.MyInt

#include("../src/Robust_evolve.jl")

nreps = 3
numinputs = 2
numoutputs = 2
numints = 6
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numints, numlevelsback=numints )
context = construct_context(numinputs)
popsize = 5
max_pop_gens = 50000
indiv_steps = 100
robust_steps = 10
fit_steps = 10
ngoals = 8   # number goals in goallist
tourn_size = 0  
csvfile = "../data/6_26/robust_evolve6_26A.csv"
gl = randgoallist(ngoals,p.numinputs,p.numoutputs)
run_robust_evolve( nreps, p, popsize, max_pop_gens, indiv_steps, robust_steps,fit_steps, gl,tourn_size, csvfile )

