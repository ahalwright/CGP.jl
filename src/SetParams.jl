
#include("../src/CGP.jl")
export p, funcs, MyInt, Ones, goallist2, goallist3
#using Base.Test

#using Main.CGP

#const MyInt = UInt16   # Also set in CGP.jl

numinputs = 2
numoutputs = 2
nodearity = 2
numinteriors = 5
numlevelsback = 8
const context = construct_contexts(numinputs)[numinputs]

goal1 = MyInt[ 0xb, 0xe ]
goal2 = MyInt[ 0x6, 0x9 ]
goal3 = MyInt[ 0xa, 0xe ]
goallist2 = [goal1, goal2,goal3]
goallist3 = [[0x50], [ 0x05], [ 0xf5], [ 0xa0], [ 0xe0], [ 0x0e], [ 0xf1], [ 0x80], [ 0x08], [ 0xfa]]

p = Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
const Ones = Main.CGP.construct_ones(numinputs)[numinputs]     


#const funcs = default_funcs()
println("SetParams: ")
print_parameters(p)

