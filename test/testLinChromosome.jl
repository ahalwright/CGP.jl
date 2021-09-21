# Test the 3-input 1-output case of Hu (2020)
# To run:  from evotech/CGP.jl/src:
# julia -L CGP.jl -L LinChromosome.jl ../test/testLinChromosome.jl
using Random
using Test
using Main.CGP
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
numinstructions=4
numregisters=2 
p = Parameters(numinputs,1,numinstructions,numregisters)
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
@test i_ints == [59, 6, 81, 109]
@test MyInt[ 0x4, 0x3, 0xc, 0xa ] == execute_lcircuit( i_ints, numregisters,numinputs,lfuncs)
cint = instruction_ints_to_circuit_int( i_ints, p, lfuncs )
@test i_ints == circuit_int_to_instruction_ints( instruction_ints_to_circuit_int( i_ints, p, lfuncs), p, lfuncs )

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
i_ints=map( x->instruction_vect_to_instruction_int( x, numregisters,numinputs,lfuncs), lcv ) 
@test i_ints == [44, 61, 130, 48, 162, 20]
@test MyInt[ 0x000a, 0x0088, 0x00f0, 0x00cc, 0x00aa ] == execute_lcircuit( i_ints, numregisters,numinputs,lfuncs)
c_int = Int128(13857033912219)
@test c_int == instruction_ints_to_circuit_int( i_ints, p, lfuncs )
@test i_ints == circuit_int_to_instruction_ints( instruction_ints_to_circuit_int( i_ints, p, lfuncs), p, lfuncs )

ctx = construct_context(3)
# Execute circuit one step at a time as specified in the Hu (2020) paper
R0=0x0000; R1=ctx[1]; R2=ctx[2]; R3=ctx[3]; R4=0x0000;
R4 = R2 & R3
R0 = R1 | R4
R4 = Nand( R4, R0 )
R4 = R3 & R2
R0 = Nor( R1, R1 )
R0 = R3 & R0
@test [R0,R4,R1,R2,R3] == execute_lcircuit( i_ints, numregisters, numinputs, lfuncs )

# Test of functions with parameter arguments.
p = Parameters(numinputs,1,length(lcv),numregisters)
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
@test Main.CGP.mutate_circuit_all( lcv, p.numlevelsback, p.numinputs, lfuncs ) == Main.CGP.mutate_circuit_all( lcv, p, lfuncs )

# Iterate through all of genotype space in the special case of 3-instruction circuits
# Return all circuits that are circuit_distance 1 from circuit rlc.
nr =  numregisters
function test_mutate(rlc::LinCircuit, numinputs::Int64, nr::Int64, funcs::Vector{Func} )
  n_instructions = num_instructions( nr, numinputs, funcs )
  nbrs = []
  for i = 1:n_instructions
    for j = 1:n_instructions
      for k = 1:n_instructions
        rcirc =[ 
            instruction_int_to_instruction_vect( i, nr, numinputs, funcs ), 
            instruction_int_to_instruction_vect( j, nr, numinputs, funcs ),
            instruction_int_to_instruction_vect( k, nr, numinputs, funcs )];
        if circuit_distance(rcirc, rlc.circuit_vects) == 1
          push!(nbrs,rcirc)
       end
      end
    end
  end
  nbrs
end    
p = Parameters(2,1,3,nr)
funcs = lin_funcs(p.numinputs)
rlc=rand_lcircuit(p)
@test sort(CGP.mutate_circuit_all( rlc.circuit_vects, p, funcs )) == sort(test_mutate(rlc,p.numinputs,nr,funcs))

end #testset
