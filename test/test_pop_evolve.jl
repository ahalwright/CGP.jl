# Test the mut_evolve_subgoal() function from src/Evolve.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = Main.CGP.MyInt

#=
# This is a test function but is included here for convenience
function test_pop_evolve(p::Parameters,popsize,funcs,max_pop_steps::Int64,max_indiv_steps::Int64,ngoals::Int64; 
    robust_sel::Bool=false)
  print_parameters(p)
  println("ngoals: ",ngoals)
  gl = randgoallist(ngoals,p.numinputs,p.numoutputs)
  (max_fit_all_gens, pop) = pop_evolve( p, popsize, max_pop_steps, max_indiv_steps, gl, funcs, robust_sel=robust_sel)
end       
=#

numinputs = 3
numoutputs = 2
popsize = 3
use_robust=true
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=10, numlevelsback=12 )
context = construct_context(numinputs)
max_indiv_steps = 20
max_pop_steps = 20
ngoals = 6   # number goals in goallist
gl = randgoallist(ngoals,p.numinputs,p.numoutputs)
robust_sel=false
all_max_sel = false
#(steps,worse,same,c,output,goallist,matched_goals,matched_goals_list)= test_pop_evolve(p,funcs,maxsteps,ngoals)
#result = test_pop_evolve(p,popsize,funcs,max_pop_steps,max_indiv_steps,ngoals,robust_sel=use_robust)
for robust_sel = true:true
  for all_max_sel = false:true
    (max_fit_all_gens, pop) = pop_evolve( p, popsize, max_pop_steps, max_indiv_steps, gl, funcs, robust_sel=robust_sel,all_max_sel=all_max_sel)
  end
end

