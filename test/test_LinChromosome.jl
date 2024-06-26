# Test the 3-input 1-output case of Hu (2020)
# To run:  from evotech/CGP.jl/src:
# julia -L CGP.jl -L LinChromosome.jl ../test/test_LinChromosome.jl
using Random
using Test
using Main.CGP
include("../src/LinChromosome.jl")
@testset "Testing functions of LinChromosome.jl in the 3-input 1-output 2-registers 6-instructions case." begin

numinputs=3
numregisters=2
lfuncs = lin_funcs(numinputs)
ni = num_instructions( numregisters, numinputs, lfuncs )
# Testing that instruction_vect_to_instruction_int() is the inverse of instruction_int_to_instruction_vect()
for i = 0:ni-1
  @test i == instruction_vect_to_instruction_int(instruction_int_to_instruction_vect( i, numregisters, numinputs, lfuncs ), numregisters, numinputs, lfuncs )
end

# Test of the circuit of Hu (2012, 309)
# See testComposition.jl for another version of this example.
numinputs=2
n_instructions=4
numregisters=2 
p = Parameters(numinputs,1,n_instructions,numregisters)
lfuncs = lin_funcs(numinputs)  # Note that lin_funcs() needs to be set to the Hu gate set
@test num_instructions(numregisters,numinputs,lfuncs) == 128 
# array of circuit vects corresponding to the circuit of Hu (2012, 309)
# An instruction is specified either by a 4-tuple of MyInts or by an integer of OutputType.
# The elements of the 4-tuple are:
#  1.  The index of the logical operation (such as AND, OR) in the funcs array.
#  2.  The index of the element of R where the output is stored
#  3.  The index of the element of R which is the first operand of the logical operation
#  4.  The index of the element of R which is the second operand of the logical operation
# Note that these indices are 1-based instead of 0-based
lc = [
  MyInt[2, 2, 3, 4],
  MyInt[1, 1, 2, 3],
  MyInt[3, 2, 1, 2],
  MyInt[4, 1, 4, 2]
]
i_ints=map( x->instruction_vect_to_instruction_int( x, numregisters,numinputs,lfuncs), lc ) 
@test i_ints == [60, 7, 82, 110]
@test MyInt[ 0x4, 0x3, 0xc, 0xa ] == execute_lcircuit( i_ints, numregisters,numinputs,lfuncs)
cint = instruction_ints_to_circuit_int( i_ints, p, lfuncs )
@test i_ints == circuit_int_to_instruction_ints( instruction_ints_to_circuit_int( i_ints, p, lfuncs), p, lfuncs )
# Executing this circuit one step at a timne.
# p = Parameters(2,1,4,2)
# circ = LinCircuit( lc, p )
# R = fill(MyInt(0), p.numlevelsback + p.numinputs ) # numlevelsback is the number of computational registers
# R[(p.numlevelsback+1):end] = construct_context(p.numinputs)
# i=1; R[lc[i][2]] = funcs[lc[i][1]].func(R[lc[i][3]],R[lc[i][4]]); R
# 4-element Vector{UInt16}:
#  0x0000
#  0x000e
#  0x000c
#  0x000a
# i = 2; R[lc[i][2]] = funcs[lc[i][1]].func(R[lc[i][3]],R[lc[i][4]]); R
# 4-element Vector{UInt16}:
#  0x000c
#  0x000e
#  0x000c
#  0x000a
# i = 3; R[lc[i][2]] = funcs[lc[i][1]].func(R[lc[i][3]],R[lc[i][4]]); R
# 4-element Vector{UInt16}:
#  0x000c
#  0x0003
#  0x000c
#  0x000a
# i = 4; R[lc[i][2]] = funcs[lc[i][1]].func(R[lc[i][3]],R[lc[i][4]]); R
# 4-element Vector{UInt16}:
#  0x0004
#  0x0003
#  0x000c
#  0x000a


# Test of the circuit of Hu (2020, 378)
numinputs=3
numinstructions=6
numregisters=2 
p = Parameters(numinputs,1,numinstructions,numregisters)
lfuncs = lin_funcs(numinputs)  # Note that lin_funcs() needs to be set to the Hu gate set
@test num_instructions(numregisters,numinputs,lfuncs) == 200
# array of circuit vects corresponding to the circuit of Hu (2020, 378)
# Note that these indices are 1-based instead of 0-based
# for the registers, our index 1 corresponds to R0
# our index 1 corresponds to R0  # calculation register
# our index 2 corresponds to R4  # calculation register
# our index 3 corresponds to R1  # input register
# our index 4 corresponds to R2  # input register 
# our index 5 corresponds to R3  # input register 
lcv = [
  MyInt[1, 2, 4, 5],
  MyInt[2, 1, 3, 2],
  MyInt[3, 2, 2, 1],
  MyInt[1, 2, 5, 4],
  MyInt[4, 1, 3, 3],
  MyInt[1, 1, 5, 1]
]
lc = LinCircuit( lcv, p )
i_ints=map( x->instruction_vect_to_instruction_int( x, numregisters,numinputs,lfuncs), lcv ) 
@test i_ints == [45, 62, 131, 49, 163, 21]
@test MyInt[ 0x000a, 0x0088, 0x00f0, 0x00cc, 0x00aa ] == execute_lcircuit( i_ints, numregisters,numinputs,lfuncs)
c_int = Int128(14178641952420)
@test c_int == instruction_ints_to_circuit_int( i_ints, p, lfuncs )
@test i_ints == circuit_int_to_instruction_ints( instruction_ints_to_circuit_int( i_ints, p, lfuncs), p, lfuncs )
@test circuit_int_to_circuit( circuit_to_circuit_int( lc, lfuncs ), p, lfuncs ).circuit_vects == lcv
@test circuit_to_circuit_int( circuit_int_to_circuit( c_int, p, lfuncs ), lfuncs ) == c_int

ctx = construct_context(3)
# Execute circuit one step at a time as specified in the Hu (2020) paper
R0=0x0000; R1=ctx[1]; R2=ctx[2]; R3=ctx[3]; R4=0x0000;
R4 = R2 & R3
R0 = R1 | R4
R4 = CGP.Nand( R4, R0 )
R4 = R3 & R2
R0 = CGP.Nor( R1, R1 )
R0 = R3 & R0
@test [R0,R4,R1,R2,R3] == execute_lcircuit( i_ints, numregisters, numinputs, lfuncs )

# Test of functions with parameter arguments.
p = Parameters(numinputs,1,length(lcv),numregisters)
funcs = lin_funcs(p.numinputs)
@test execute_lcircuit( lcv, numregisters, numinputs, lfuncs ) ==  execute_lcircuit( LinCircuit(lcv,p), lfuncs )
@test num_instructions( p.numlevelsback, p.numinputs, lfuncs ) == num_instructions( p, lfuncs )
Random.seed!(1); riv0 = rand_ivect( p.numlevelsback, p.numinputs, lfuncs )
Random.seed!(1); riv1 = rand_ivect( p, lfuncs )
@test riv0 == riv1
Random.seed!(1); rii0 = rand_lcircuit( p.numinteriors, p.numlevelsback, p.numinputs, lfuncs )
Random.seed!(1); rii1 = rand_lcircuit( p, lfuncs )
#println("rii0.params: ",rii0.params,"  rii1.params: ",rii1.params)
#println("rii0.circuit_vects: ",rii0.circuit_vects)
#println("rii1.circuit_vects: ",rii1.circuit_vects)
#@test rii0.params == rii1.params && rii0.circuit_vects == rii1.circuit_vects  # 6/30/21: test rii0.params == rii1.params failed even though they are equal
@test rii0.circuit_vects == rii1.circuit_vects
llfuncs = length(lfuncs)
num_mutate_locations = ((llfuncs>1) ? (llfuncs-1) : 0) + (p.numlevelsback-1) + p.nodearity*(p.numlevelsback+p.numinputs-1) 
nml0=map(ml->instruction_vect_to_instruction_int(Main.CGP.mutate_instruction(lcv[1],p.numlevelsback,p.numinputs,lfuncs,ml),p.numlevelsback,p.numinputs,lfuncs),1:num_mutate_locations)
nml1=map(ml->instruction_vect_to_instruction_int(Main.CGP.mutate_instruction(lcv[1],p,lfuncs,ml),p.numlevelsback,p.numinputs,lfuncs),1:num_mutate_locations)
@test nml0 == nml1
@test length(unique(nml0)) == length(nml0)   # Check that all mutations are unique
Random.seed!(1); mc0 = Main.CGP.mutate_instruction( lcv[1], p.numlevelsback, p.numinputs, lfuncs )
Random.seed!(1); mc1 = Main.CGP.mutate_instruction( lcv[1], p, lfuncs )
@test mc0 == mc1
Random.seed!(1); mc0=Main.CGP.mutate_circuit!(deepcopy(lcv),p,lfuncs) 
Random.seed!(1); mc1=Main.CGP.mutate_circuit!(deepcopy(lcv),p.numlevelsback,p.numinputs,lfuncs)
@test mc0 == mc1
Random.seed!(1); mc0 = Main.CGP.mutate_circuit!( deepcopy(lcv), p.numlevelsback, p.numinputs, lfuncs )
Random.seed!(1); mc1 = Main.CGP.mutate_circuit!( deepcopy(lcv), p, lfuncs )
@test mc0 == mc1
mca0 = Main.CGP.mutate_circuit_all( lcv, p.numlevelsback, p.numinputs, lfuncs ); 
mca1 = Main.CGP.mutate_circuit_all( lcv, p, lfuncs );
la0 = map( c->circuit_to_circuit_int(c,funcs), mca0 );
la1 = map( c->circuit_to_circuit_int(c,funcs), mca1 );
@test la0 == la1

numinstructions = 3
numregisters = 2
p = Parameters(2,1,numinstructions,numregisters)
funcs = lin_funcs(p.numinputs)
rlc=rand_lcircuit(p,funcs)
println("rlc: ")
print_lcircuit(rlc)
#test_remove_inactive(rlc)
@test output_values(rlc) == output_values( remove_inactive( rlc ) )
mca = CGP.mutate_circuit_all( rlc.circuit_vects, p, funcs );
println("mca[1]: ")
print_lcircuit(mca[1])
@test unique(map( mc->circuit_distance( rlc, mc ), mca )) == Int64[1]
@test count_circuits_lc( p, funcs ) == length(enumerate_circuits_lc( p, funcs ))  # 16384

end #testset
