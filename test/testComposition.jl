using Test
using Main.CGP

@testset "compositional" begin  # @testset makes random_chromosome(p) always return same chromosome

p = Parameters( 2, 1, 3, 3 )
funcs = default_funcs(p.numinputs)
ctx = construct_context(p.numinputs)

cg3 = CGP.Cg( NAND, [1,2], 1 )
cg4 = CGP.Cg( XOR, [2,3], 1 )
cg5 = CGP.Cg( AND, [4,3], 1 )
cg6 = CGP.Cg( NOR, [3,5], 1 )
cg7 = CGP.Cg( NOR, [5,4], 1 )
cc7 = CGP.CC( CGP.CC[cg3, cg4, cg5, cg6, cg7], [1,2], 1)
ch7 = circuit((1,2), ((3,NAND,1,2), (4,XOR,2,3), (5,AND,4,3), (6,NOR,3,5), (7,NOR,5,4)))
@test execute!(cg3,ctx,active_only=true) == execute!(cg3,ctx,active_only=false) == 0x0007
@test execute!(cg4,[0x000c, 0x000a, 0x0007],active_only=true) == execute!(cg4,[0x000c, 0x000a, 0x0007],active_only=false) == 0x000d
@test execute!(cg5,[0x000c, 0x000a, 0x0007, 0x000d],active_only=true) == execute!(cg5,[0x000c, 0x000a, 0x0007, 0x000d],active_only=false) == 0x0005
@test execute!(cg6,[0x000c, 0x000a, 0x0007, 0x000d, 0x0005],active_only=true) == execute!(cg6,[0x000c, 0x000a, 0x0007, 0x000d, 0x0005],active_only=false) == 0x0008
@test execute!(cg7,[0x000c, 0x000a, 0x0007, 0x000d, 0x0005, 0x0008],active_only=true) == execute!(cg7,[0x000c, 0x000a, 0x0007, 0x000d, 0x0005, 0x0008],active_only=false) == 0x0002
@test output_values(ch7)[1] == execute!(cc7,ctx) == execute!(cg7,[0x000c, 0x000a, 0x0007, 0x000d, 0x0005, 0x0008],active_only=true)

cc1 = CGP.CC( CGP.CC[cg3, cg4], [1,2], 1 )
@test execute!(cc1,ctx,active_only=true) == execute!(cc1,ctx,active_only=false) == 0x000d
cc4 = CGP.CC( CGP.CC[cg4, cc1], [2,1], 1 )
@test execute!(cc4,ctx,active_only=true) == execute!(cc4,ctx,active_only=false) == 0x000b
cc6 = CGP.CC( CGP.CC[cc1,cc4], [1,2], 1 )
@test execute!(cc6,ctx,active_only=true) == execute!(cc6,ctx,active_only=false) == 0x000b

numinputs = 3
p = Parameters(numinputs,1,4,3)
ctx = construct_context(numinputs)
ccg = CGP.Cg(AND,[2,3],1)
@test execute!(ccg,ctx) == 0x0088
cc2 = CGP.CC(CGP.CC[CGP.Cg(AND,[2,3],1),CGP.Cg(XOR,[1,4],1)],[1,2,3],1)
@test execute!(cc2,ctx) == 0x0078

ch = random_chromosome(p)
print_circuit(ch)
println("output_values(ch)[1]: ",@sprintf("0x00%x",output_values(ch)[1]))
println("chromosome_to_circuit(ch): ",chromosome_to_circuit(ch))
@test execute!(chromosome_to_circuit(ch),ctx) == output_values(ch)[1]

cg4 = CGP.Cg( NAND, [1,2], 1 ) #4
cg5 = CGP.Cg( XOR, [4,3], 1 )  #5
cg6 = CGP.Cg( AND, [4,3], 1 )  #6
cg7 = CGP.Cg( NOR, [3,5], 1 )  #7
cg8 = CGP.Cg( OR, [5,6] ,1 )  #8
cc1 = CGP.CC( CGP.CC[cg4, cg5 ], [2,1,3], 1 )
cc2 = CGP.CC( CGP.CC[cg5, cg4, cg7], [2,3,1], 1 )
cc3 = CGP.CC( CGP.CC[cg4, cg5, cg6, cg7, cg8], [3,1,2], 1 )
cc4 = CGP.CC( CGP.CC[cc1], [1,3,2], 1 )
cc5 = CGP.CC( CGP.CC[cc4, cc2], [1,3,2], 1 )  
@test execute!(cg4,ctx) == 0x003f
@test execute!(cc1,ctx) == 0x0095 == execute_chromosome(circuit_to_chromosome(cc1),ctx)[1] 
@test execute!(cc2,ctx) == 0x0008 == execute_chromosome(circuit_to_chromosome(cc2),ctx)[1] 
@test execute!(cc3,ctx) == 0x00df == execute_chromosome(circuit_to_chromosome(cc3),ctx)[1]
@test execute!(cc4,ctx) == 0x0093 == execute!(cc1,ctx[[1,3,2]]) 
@test execute!(cc5,ctx) == 0x0008

@test execute_chromosome(ch,ctx)[1] == execute!(chromosome_to_circuit(ch),ctx)
@test execute_chromosome(circuit_to_chromosome(cc1),ctx)[1] == execute!(cc1,ctx)
cc9 = random_gate_circuit( 3, 4, funcs)
ch9 = circuit_to_chromosome(cc9); print_circuit(ch9)
@test execute_chromosome(ch9,ctx)[1] == execute!(chromosome_to_circuit(ch9),ctx)
@test execute_chromosome(circuit_to_chromosome(cc9),ctx,permute_context=true)[1] == execute!(cc9,ctx)
@test execute_chromosome(circuit_to_chromosome(cc9),ctx)[1] == execute!(cc9,ctx)

Random.seed!(1);
needs_list = [([0x000a],2)]
numinputs=2
numcompositions=2
numcircuits=4
funcs = [AND,OR,XOR,NAND,NOR,IN1,IN2]
(needs,ac) =comp_experiment( numinputs, numcompositions, needs_list, numcircuits, funcs );
@test needs[1].goal[1] == 0x000a  # Might not work if random number generation changes

end  # @testset

