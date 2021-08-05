using Test
numinputs = 3
p = Parameters(numinputs,1,4,3)
ch = random_chromosome(p)

@testset "compositional" begin  # @testset makes random_chromosome(p) always return same chromosome
ctx = construct_context(numinputs)
ccg = CompGate([2,3],AND,1)
@test execute(ccg,ctx) == 0x0088
cc2 = CompCircuit([1, 2, 3], [CompGate([2,3],AND,1),CompGate([1,4],XOR,1)],1)
@test execute(cc2,ctx) == 0x0078

print_circuit(ch)
println("output_values(ch)[1]: ",@sprintf("0x00%x",output_values(ch)[1]))
println("comp_circuit(ch): ",comp_circuit(ch))
@test execute(comp_circuit(ch),ctx) == output_values(ch)[1]
end
