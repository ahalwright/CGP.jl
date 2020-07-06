# Test the mut_evolve_subgoal() function from src/Evolve.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = Main.CGP.MyInt

# This is a test function but is included here for convenience
function test_mut_evolve(p::Parameters,funcs,max_steps::Int64,ngoals::Int64; use_robustness::Bool=false)
  #println("ngoals: ",ngoals)
  c = random_chromosome(p,funcs)
  gl = randgoallist(ngoals,c.params.numinputs,c.params.numoutputs)
  (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, gl, funcs, max_steps,  use_robustness=use_robustness)
  @assert output == output_values(c)
  (c,step,worse,same,better,output,gl,matched_goals,matched_goals_list)
end       

numinputs = 2
numoutputs = 2
use_robust=false
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=10, numlevelsback=12 )
context = construct_context(numinputs)
c = random_chromosome( p, funcs )
maxsteps = 50000
ngoals = 2   # number goals in goallist
#(steps,worse,same,c,output,goallist,matched_goals,matched_goals,matched_goals_list)= test_mut_evolve(p,funcs,maxsteps,ngoals)
(c,steps,worse,same,better,output,gl,matched_goals,matched_goals_list)= test_mut_evolve(p,funcs,maxsteps,ngoals,use_robustness=use_robust)
#println("result: ",(c,steps,worse,same,better,output,gl,matched_goals_list))

