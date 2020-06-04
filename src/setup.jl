include("CGP.jl")
#include("Avg_mut_robustness.jl")
#include("nextc")
#include("Subgoal_evolution.jl")
numinputs = 2
funcs = default_funcs(numinputs)
p = Main.CGP.Parameters( numinputs=numinputs, numoutputs=2, numinteriors=5, numlevelsback=8 )
context = construct_context(numinputs)
c = random_chromosome(p,funcs)
goallist = randgoallist(5,p.numinputs,p.numoutputs)
MyInt = Main.CGP.MyInt

