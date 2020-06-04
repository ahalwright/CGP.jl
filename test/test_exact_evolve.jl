# Test the mut_evolve_exact() function from src/Evolve.jl

#include("../src/CGP.jl")
using Main.CGP
const MyInt = UInt8  

goal1 = Main.CGP.MyInt[ 0xb, 0xe ]
goal2 = Main.CGP.MyInt[ 0x6, 0x9 ]
goal3 = Main.CGP.MyInt[ 0xa, 0xe ]
#goallist = [goal1, goal2,goal3]

goallist = [[0x0,0x1],[0x07, 0x01],[0x0d, 0x0e],[0x07, 0x0d],[0xb,0xe],[0x6,0x9],[0xa,0xe]]

numinputs = 2
funcs = default_funcs(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=2, numinteriors=6, numlevelsback=6 )
context = construct_context(numinputs)
c = random_chromosome( p, funcs )
max_steps = 20
result = mut_evolve_exact( random_chromosome( p, funcs ), context, goallist, max_steps )

Random.seed!(5);
c = random_chromosome( p, funcs )
ec = eval_chromosome(ch5,goallist)
println("ec: ",ec)
@assert ec[1] == [0, 0, 2, 0, 1, 0, 1]
@assert chrome_check(c,goallist,0) == 1
@assert chrome_check(c,goallist,1) == 1
@assert chrome_check(c,goallist,2) == 0
@assert chrome_check(c,goallist,3) == -1
