# Test the mut_evolve_subgoal() function from src/Evolve.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = UInt8
#include("../src/Evolve.jl")   # included in CGP.jl

# This is a test function but is included here for convenience
function mut_evolve(p,funcs,max_steps,goallist_length)
  println("goallist_length: ",goallist_length)
  c = random_chromosome(p,funcs)
  gl = randgoallist(goallist_length,c.params.numoutputs)
  (step,worse,same,c,output,goallist,matched_goals,matched_goals_list) = mut_evolve_subgoal( c, gl, funcs, max_steps )
  @assert output == output_values(c)
  ((step,worse,same),c,gl,output,goallist,matched_goals,matched_goals_list)
end       

numinputs = 2
numoutputs = 1
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=4, numlevelsback=12 )
context = construct_context(numinputs)
c = random_chromosome( p, funcs )
nsteps = 40000
ngoals = 2   # number goals in goallist
#result = mut_evolve_subgoal(p,funcs,nsteps,ngoals)
#println("length result: ",length(result),"  result: ",result)
#(steps,worse,same,c,output,goallist,matched_goals,matched_goals_list)= mut_evolve(p,funcs,nsteps,ngoals)
((steps,worse,same),c,output,goallist,matched_goals,matched_goals_list)= mut_evolve(p,funcs,nsteps,ngoals)
println("result: ",((steps,worse,same),c,output,goallist,matched_goals,matched_goals_list))

#=
goal1 = Main.CGP.MyInt[ 0xb, 0xe ]
goal2 = Main.CGP.MyInt[ 0x6, 0x9 ]
goal3 = Main.CGP.MyInt[ 0xa, 0xe ]
#goallist = [goal1, goal2,goal3]
goallist = [[0x0,0x1],[0x07, 0x01],[0x01, 0x0d],[0x07, 0x0d]] 
outputs =  [0x07, 0x0d] 
=#
