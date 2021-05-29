# Test the mut_evolve() function from src/Evolve.jl
# Worked 5/29/21

#include("../src/CGP.jl")
using Main.CGP
MyInt = Main.CGP.MyInt

# This is a test function but is included here for convenience
function test_mut_evolve(p::Parameters, funcs::Vector{Func}, max_steps::Int64, goallist::GoalList; use_robustness::Bool=false,
      avgfitness::Bool=false )
  #println("ngoals: ",ngoals)
  c = random_chromosome(p,funcs)
  println("gl: ",gl)
  fit_limit = Float64(c.params.numoutputs)
  (c,step,worse,same,better,output,matched_goals,matched_goals_list) = 
      mut_evolve( c, gl, funcs, max_steps,  use_robustness=use_robustness, avgfitness=avgfitness, fit_limit=fit_limit, print_improvements=false ) 
  if step < max_steps
    @assert output == output_values(c)
  end
  (c,step,worse,same,better,output,gl,matched_goals,matched_goals_list)
end       

avgfitness = true
numinputs = 4
numoutputs = 1
numinteriors = 10
use_robust=false
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numinteriors, numlevelsback=numinteriors )
context = construct_context(numinputs)
c = random_chromosome( p, funcs )
print_circuit(c)
maxsteps = 200000
ngoals = 1   # number goals in goallist
gl = randgoallist(ngoals,c.params.numinputs,c.params.numoutputs)
#(steps,worse,same,c,output,goallist,matched_goals,matched_goals,matched_goals_list) = 
#      test_mut_evolve(p,funcs,maxsteps,ngoals, avgfitness=avgfitness)
(ec,steps,worse,same,better,output,gl,matched_goals,matched_goals_list ) = 
      test_mut_evolve(p,funcs,maxsteps,gl,use_robustness=use_robust, avgfitness=avgfitness )
print_circuit(ec)
#(nc,steps) = neutral_evolution(c,gl[1],maxsteps)
#print_circuit(nc)
#println("result: ",(c,steps,worse,same,better,output,gl,matched_goals_list))
#(c,steps,worse,same,better,output,gl,matched_goals_list)

