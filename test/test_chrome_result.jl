# Test Chrome_result.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = UInt8
include("../src/Chrome_result.jl")

randgoal(numoutputs) = [rand(collect(0x0:0xf)) for _ in 1:numoutputs]  # Assumes MyInt = UInt8
randgoallist(len,numouputs) = [randgoal(numoutputs) for _ = 1:len]

function mut_evolve(p,funcs) 
  c = random_chromosome(p,funcs)
  gl = randgoallist(6,3)
  (worse,same,c,output,matched_goals,matched_goals_list) = mut_evolve_subgoal( c, gl, funcs, 500 )
  (c,gl,output,matched_goals,matched_goals_list)
  @assert output == output_values(c)
end
#goal1 = Main.CGP.MyInt[ 0xb, 0xe ]
#goal2 = Main.CGP.MyInt[ 0x6, 0x9 ]
#goal3 = Main.CGP.MyInt[ 0xa, 0xe ]
#goallist = [goal1, goal2,goal3]

#goallist = [[0x0,0x1],[0x07, 0x01],[0x01, 0x0d],[0x07, 0x0d]] 
#output =  [0x07, 0x0d] 

numcomponents = 6
routput = rand([0x0,0x1,0x2,0x3],numcomponents)
println("routput: ",routput)
rgoal =   rand([0x0,0x1,0x2,0x3],numcomponents)
println("rgoal: ",rgoal)
cm = components_matched( routput, rgoal )
println("number_components_matched: ",length(cm))
println("components_matched: ",cm)
check = map(x->(routput[cm[x][2]],rgoal[cm[x][3]]), collect(1:length(cm)))
println("check: ",check)
println(map(x->(x[1]==x[2]),check),"  all components should be 1")

numgoals = 4
rgoallist = [rand([0x0,0x1,0x2,0x3],numcomponents) for _=1:numgoals]
gm = goals_matched(routput,rgoallist)
println("goals matched: ",gm)



#=
numinputs = 2
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=2, numinteriors=6, numlevelsback=6 )
context = construct_context(numinputs)
c = random_chromosome( p, funcs )
max_steps = 20
#result = mut_evolve_exact( random_chromosome( p, funcs ), context, goallist, max_steps )
=#
