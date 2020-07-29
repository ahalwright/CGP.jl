# Simple test of fault tolerance functions

using Test
using Main.CGP
MyFunc = Main.CGP.MyFunc
funcs = default_funcs(3)   
Not = Main.CGP.Not  # Not() only inverts the relevant bits
c = build_chromosome(
  (1,2,3),
  ((NAND,Integer[3, 1]),(NOR,Integer[2, 1]),(OR,Integer[4, 5])),
  (5,6))

context = construct_context(c.params.numinputs)
outv = output_values(c)

@testset "test correctness of execute_chromomsome_ft() and fault_tolerance_fitness()" begin
# Compute output for c
n1 = context[1]
n2 = context[2]
n3 = context[3]
n4 = Nand(n3,n1)
n5 = Nor(n2,n1)
n6 = n4 | n5
n7 = n5
n8 = n6
@test outv == [n7,n8]

# compute output for c with gate 6 negated
p6 = Not(n6)
p7 = n5
p8 = p6
e6 = execute_chromosome_ft( c, context, 6 );
@test  e6 == [p7,p8]

# compute output for c with gate 5 negated
q5 = Not(n5)
q6 = n4 | q5
q7 = q5
q8 = q6
e5 = execute_chromosome_ft( c, context, 5 );
@test  e5 == [q7,q8]

# compute output for c with gate 4 negated
r4 = Not(n4)
r5 = Nor(n2,n1)
r6 = r4 | r5
r7 = r5
r8 = r6
e4 = execute_chromosome_ft( c, context, 4 );
@test  e4 == [r7,r8]

# compute output for c with gate 3 negated
s3 = Not(n3)
s4 = Nand(s3,n1)
s5 = Nor(n2,n1)
s6 = s4 | s5
s7 = s5
s8 = s6
e3 = execute_chromosome_ft( c, context, 3 );
@test  e3 == [s7,s8]

funcs = default_funcs(2)  # Reset for chromosome c1 which has 2 inputs
# Flipping the output of gate 3 of c1 changes the output from 0x00 to 0x0f
c1 = build_chromosome( (1,2), ((ZERO,Integer[]),), (3,))
println("fault_tolerance_fitness(c1): ",fault_tolerance_fitness(c1))
@test isapprox(fault_tolerance_fitness(c1),0.0)

funcs = default_funcs(3)  # Reset for chromosome nulc which has 3 inputs
# Flipping the output of gate 4 (the first interior node) has no effect on the output
# Flipping the output of gate 5 (the second interior node)  flips the first output but not the second
# Flipping the output of gate 6 (the third interior node)  flips the second output but not the first
nulc = build_chromosome(
  (1,2,3),
  ((NAND,Integer[3, 1]),(ZERO,Integer[]),(ONE,Integer[])),
  (5,6))

@test isapprox(fault_tolerance_fitness(nulc),2/3)
end  # @testset
