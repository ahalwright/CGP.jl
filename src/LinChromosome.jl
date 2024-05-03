# Linear GP circuit chromosome functions
# When Parameters are used, numinteriors specifies numinstructions.
# When Parameters are used, numlevelsback specifies numregisters.
# A circuit is a sequence of instructions on the register array R 
#   which has length numregisters + numinputs
# R is normally local to the function where is is used.
# The first numregisters components of R are calculation registers.
# The remaining numinputs components are read-only input registers which are 
#   initialized with the components of the current context.
# Even though there is a nodearity keyword value, this file assumes nodearity==2
# An instruction is specified either by a 4-tuple of MyInts or by an integer of OutputType.
# The elements of the (2+nodearity)-tuple are (assuming nodearity==2):
#  1.  The index of the logical operation (such as AND, OR) in the funcs array.
#  2.  The index of the element of R where the output is stored
#  3.  The index of the element of R which is the first operand of the logical operation
#  4.  The index of the element of R which is the second operand of the logical operation
# The functions vect_to_int() and int_to_vect() convert between these representations.
# When parameters are used, numinputs and numoutputs have the same meaning for the cartesian and linear representatons.
# numinteriors means numinstructions for the linear representation
# numlevelsback means numregisters for the linear representation
# There are 3 alternative representations for the circuit_vects component of a  LinCicuit.
#   1.  A list of numinstructions vectors where each vector is  (2+nodearity)-tuple described above.
#   2.  A list of numinstructions Int64s where each of these is a 1-based encoding of the vectors
#   3.  A Int128 integer which is an encoding of the integers of step 3.
# See notes/3_24.txt for outline
# See test/testLinCircuit.jl for tests.
#=
export LinCircuit, output_values
export execute_lcircuit, instruction_vect_to_instruction_int, instruction_int_to_instruction_vect, lincomplexity, print_lcircuit 
export vect_to_int, int_to_vect, rand_ivect, execute_random_circuit, degeneracy, recover_phenotype
export num_instructions # This is the total number of possible instructions for a parameter setting
export rand_lcircuit, mutate_circuit!, mutate_instruction, mutate_circuit_all, mutate_all, print_circuit
export instruction_vects_to_instruction_ints, instruction_ints_to_instruction_vects
export instruction_ints_to_circuit_int, circuit_int_to_instruction_ints
export circuit_int_to_circuit, circuit_to_circuit_int, enumerate_circuits_lc, count_circuits_lc
export count_genotypes_lc, count_genotypes_lc_mt
export number_active, remove_inactive
=#
OutputType = Int64
LinCircuit = CGP.LinCircuit

#=
mutable struct LinCircuit
    # Here, numinteriors represents number of gates, numlevelsback represents number of computational registers
    # the First numoutputs registers are the output registers
    circuit_vects::Vector{Vector{MyInt}}
    params::Parameters   
end
=#

# Print an LinCircuit in the format used in Hu and Banzhaf papers.
function print_lcircuit( lc::LinCircuit, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = lin_funcs( lc.params.numinputs )
  end
  println("funcs: ",funcs)
  p = lc.params
  c = lc.circuit_vects
  #for i = 1:p.numinteriors  # p.numinteriors is the number of instructions
  for i = 1:length(lc.circuit_vects)  
    println( "R$(c[i][2]) = R$(c[i][3]) $(funcs[[c[i][1]]][1].name) R$(c[i][4])")
  end
end

# Returns the phenotype (Goal) computed by c
function output_values( c::LinCircuit, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs( c.params.numinputs )
  end
  #println("funcs: ",funcs)
  R = execute_lcircuit( c, funcs )
  R[1:c.params.numoutputs]
end  

function output_values( cv::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func} )
  output_values( LinCircuit( cv, p ), funcs )
end

# ci is a circuit_int
function output_values( ci::Integer, p::Parameters, funcs::Vector{Func} )
  output_values(instruction_ints_to_instruction_vects(circuit_int_to_instruction_ints(ci,p,funcs),p,funcs),p,funcs)
end

# Not used.  Assumes that R and funcs are in the execution environment
function vect_to_funct( lc::Vector{MyInt} )
  (lc)->R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]]) 
end

# Executes the circuit, returns the register array R.
# numregisters is the number of calculation registers, not the total number of registers
# The outputs are  R[1:numoutputs]
function execute_lcircuit( circuit_ints::Vector{OutputType}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  lcv = map( ci->instruction_int_to_instruction_vect( ci, numregisters, numinputs, funcs ), circuit_ints )
  execute_lcircuit( lcv, numregisters, numinputs, funcs )
end

function execute_lcircuit( circuit::LinCircuit, funcs::Vector{Func} )
  p = circuit.params
  R = fill(MyInt(0), p.numlevelsback + p.numinputs ) # numlevelsback is the number of computational registers
  R[(p.numlevelsback+1):end] = construct_context(p.numinputs)  
  for lc in circuit.circuit_vects
    #println("lc: ",lc,"  R: ",R)
    R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
  end
  R
end

# Executes the circuit, returns the register array R.
# The outputs are  R[1:numoutputs]
function execute_lcircuit( circuit_vects::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  R = fill(MyInt(0), numregisters+numinputs )
  R[numregisters+1:end] = construct_context(numinputs)
  for lc in circuit_vects
    #println("lc: ",lc,"  func: ",funcs[lc[1]],"  R: ",R)
    R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
    #println("lc: ",lc,"  func: ",funcs[lc[1]],"  R: ",R)
  end
  R
end

# Returns the Tononi complexity of the circuit
function lincomplexity( circuit_vects::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  p = Parameters( numinputs, 1, length(circuit_vects), numregisters )
  lincomplexity( LinCircuit(circuit_vects,p), funcs, nodearity=nodearity )
end

# Returns the Tononi complexity of the circuit
function lincomplexity( circuit::LinCircuit, funcs::Vector{Func}; nodearity::Int64=2 )
  p = circuit.params
  R = fill(MyInt(0), p.numlevelsback+p.numinputs )
  X = fill(MyInt(0), length(circuit.circuit_vects) )  
  R[p.numlevelsback+1:end] = construct_context(p.numinputs)
  i = 1
  for lc in circuit.circuit_vects
    #println("lc: ",lc,"  func: ",funcs[lc[1]],"  R: ",R)
    #println("lc: ",lc,"  R: ",R)
    R[lc[2]] = X[i] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
    #println("lc: ",lc,"  func: ",funcs[lc[1]],"  R: ",R)
    i += 1
  end
  complexity5( X, p.numinputs )
end

# The number of possible instructions with these settings.
# This is distinct from numinstructions which is the number of instructions in a circuit.
function num_instructions( numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; 
     nodearity::Int64=2 )
  length(funcs)*numregisters*(numinputs+numregisters)^nodearity
end

function num_instructions( p::Parameters, funcs::Vector{Func} )
  num_instructions( p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

#  Assumes 1 output register 
#  The calculation registers are 1 and 2 and the output register is 1.
#  Registers 3, 4, 5 are the input registers.
#  Instruction fields are 1 for function, 2 for output, 3 and 4 for input
function number_active( circ::LinCircuit; return_inst_active::Bool=false )
  reg_active = Bool[true, false, false, false, false]  # Initially only register field 1 is active.
  inst_active = fill(false,circ.params.numinteriors)   # circ.params.numinteriors is the number of instrctions
  for i = circ.params.numinteriors:-1:1   # Iterate through instructions in reverse order
    inst = circ.circuit_vects[i]
    if reg_active[inst[2]] 
      inst_active[i] = true
      if inst[3] in [MyInt(1),MyInt(2)]
        reg_active[inst[3]] = true
      end
      if inst[4] in [MyInt(1),MyInt(2)]
        reg_active[inst[4]] = true
      end
      if !reg_active[inst[3]] && !reg_active[inst[4]]
        reg_active[inst[2]] = false
      end
    end
  end
  if !return_inst_active
    sum(inst_active)
  else
    (sum(inst_active), reg_active, inst_active )
  end
end

# Removes the inactive instructions from circ
function remove_inactive( circ::LinCircuit, active::Vector{Bool} )
  active_indices = filter( x->x!=0, [ active[i] ? i : 0 for i = 1:length(active) ] )
  new_params = Parameters( circ.params.numinputs, circ.params.numoutputs, length(active_indices), circ.params.numlevelsback )
  LinCircuit( circ.circuit_vects[active_indices], new_params )
end

# Removes the inactive instructions from circ
function remove_inactive( circ::LinCircuit )
  (num_active, reg_active, active ) = number_active( circ, return_inst_active=true )
  remove_inactive( circ, active )
end

# Returns nothing if successful
function test_remove_inactive( circ::LinCircuit )
  @assert output_values(circ) == output_values( remove_inactive( circ ) )
end

function validate_lcircuit( circ::LinCircuit )
  @assert length(circ.circuit_vects) == circ.params.numinteriors
end
  
# Convert a vectors of integers to a unique integer.
# The components of the vector are 1-based rather than 0-based.
# Example 1:  
# julia> vect_to_int( [3,4,5,6],[10,10,10,10])
#  2345
# Example 2:
# julia> vect_to_int([3,2,1,4],[7,2,5,8])
#  203
function vect_to_int( ints::Vector{Int64}, max::Vector{Int64} )
  result = Int128(ints[1])-1
  for i = 2:length(ints)
    result *= max[i]
    result += ints[i]-1
  end
  result
end

# julia> int_to_vect( 2345, [10,10,10,10])
# 4-element Array{Int64,1}:
#  Int64[ 3, 4, 5, 6 ]
# julia> int_to_vect(203,[7,2,5,8])
# 4-element Array{Int64,1}:
#  Int64[ 3, 2, 1, 4] 
function int_to_vect( int::Integer, max::Vector{Int64} )
  ints = zeros(Int64,length(max))
  for i = length(max):-1:1
    ints[i] = int % max[i] + 1
    int = div(int,max[i])
  end
  ints
end
    
# Converts an instruction specified by a vector to the instruction specified by an integer.
function instruction_vect_to_instruction_int( inst_vect::Vector{MyInt}, numregisters::Int64, numinputs::Int64, 
    funcs::Vector{Func}; nodearity::Int64=2 )
  result = Int64(inst_vect[1]-1)  # Function code
  multiplier = 1
  #println("multiplier: ",multiplier,"  result: ",result)
  multiplier = numregisters
  result *= multiplier
  result += inst_vect[2]-1
  #println("multiplier: ",multiplier,"  result: ",result)
  multiplier = numregisters+numinputs
  for j = 1:nodearity
    result *= multiplier
    result += inst_vect[2+j]-1
    #println("multiplier: ",multiplier,"  result: ",result)
  end
  result + 1    # Convert 1-based from 0-based
end

# Converts an instruction specified by an integer to the instruction specified by a vector.
function instruction_int_to_instruction_vect( inst_int::OutputType, numregisters::Int64, numinputs::Int64, 
    funcs::Vector{Func}; nodearity::Int64=2 )  
  #println("ii_to_iv inst_int: ",inst_int)
  result = fill(MyInt(0),2+nodearity)
  multiplier = numregisters+numinputs
  inst_int -= 1   # Convert from 1-based to 0-based so that mod and div work
  for j = nodearity:-1:1
    try
      result[j+2] = inst_int % multiplier + 1
    catch
      println("catch inst_int: ",inst_int,"  multiplier: ",multiplier)
      error(" error ")
    end
    inst_int = div(inst_int,multiplier)
    #println("j+2: ",j+2,"  multiplier: ", multiplier,"  inst_int: ",inst_int,"  result[j+2]: ",result[j+2])
  end
  multiplier = numregisters
  result[2] = inst_int % multiplier + 1
  inst_int = div(inst_int,multiplier)
  #println("i: ",2,"  multiplier: ", multiplier,"  inst_int: ",inst_int,"  result[2]: ",result[2])
  multiplier = length(funcs)
  result[1] = inst_int % multiplier + 1
  inst_int = div(inst_int,multiplier)
  #println("i: ",1,"  multiplier: ", multiplier,"  inst_int: ",inst_int,"  result[1]: ",result[1])
  result
end 

function instruction_vects_to_instruction_ints( cv::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; 
    nodearity::Int64=2 )
  map(x->instruction_vect_to_instruction_int(x, numregisters, numinputs, funcs, nodearity=nodearity ), cv )
end

function instruction_vects_to_instruction_ints( circ::LinCircuit, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs( circ.params.numinputs )
  end
  instruction_vects_to_instruction_ints( circ.circuit_vects, circ.params.numlevelsback, circ.params.numinputs, funcs )
end

function instruction_ints_to_instruction_vects( c_ints::Vector{Int64}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func};
    nodearity::Int64=2 )
  map(x->instruction_int_to_instruction_vect(x, numregisters, numinputs, funcs, nodearity=nodearity ), c_ints )
end

function instruction_ints_to_instruction_vects( c_ints::Vector{Int64}, p::Parameters, funcs::Vector{Func}=Func[];
    nodearity::Int64=2 )
  if length(funcs) == 0
    funcs = default_funcs( p.numinputs )
  end
  instruction_ints_to_instruction_vects( c_ints, p.numlevelsback, p.numinputs, funcs )
end

#=
function circuit_vect_to_circuit_ints( c::LinCircuit, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs( c.params.numinputs )
  end
  circuit_vect_to_circuit_ints( c.circuit_vects, c.params.numinteriors, c.params.numinputs, funcs )
end
=#

function instruction_ints_to_circuit_int( c_ints::Vector{Int64}, p::Parameters, funcs::Vector{Func} )
  numregisters = p.numlevelsback
  multiplier = num_instructions( p, funcs )
  #println("multiplier: ",multiplier)
  vect_to_int( c_ints, fill(multiplier,length(c_ints)))
end

function circuit_int_to_instruction_ints( c_int::Integer, p::Parameters, funcs::Vector{Func} ) 
  numinstructions = p.numinteriors
  multiplier = num_instructions( p, funcs )
  #println("multiplier: ",multiplier)
  int_to_vect( c_int, fill( multiplier, numinstructions ))
end

function circuit_int_to_instruction_ints( c_int::Integer, p::Parameters, funcs::Vector{Func} ) 
  numinstructions = p.numinteriors
  multiplier = num_instructions( p, funcs )
  #println("multiplier: ",multiplier)
  int_to_vect( c_int, fill( multiplier, numinstructions ))
end

function circuit_to_circuit_int( circ::LinCircuit, funcs::Vector{Func} )
  instruction_ints_to_circuit_int( instruction_vects_to_instruction_ints( circ, funcs ), circ.params, funcs )
end

function Lincircuit_int_to_circuit( c_int::Integer, p::Parameters, funcs::Vector{Func} )
  i_vects = instruction_ints_to_instruction_vects( circuit_int_to_instruction_ints( c_int, p, funcs), p, funcs )
  LinCircuit( i_vects, p )
end

function circuit_int_to_circuit( c_int::Integer, p::Parameters, funcs::Vector{Func} )
  i_vects = instruction_ints_to_instruction_vects( circuit_int_to_instruction_ints( c_int, p, funcs), p, funcs )
  LinCircuit( i_vects, p )
end

# Random instruction  
# Assumes gates are arity 2
function rand_ivect( numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  MyInt[ rand(1:length(funcs)), rand(1:numregisters),  rand(1:(numregisters+numinputs)), rand(1:(numregisters+numinputs)) ]
end

# Random instruction  
# Assumes gates are arity 2 which means that instructions have 2 components
function rand_ivect( p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )
  MyInt[ rand(1:length(funcs)), rand(1:p.numlevelsback),  rand(1:(p.numlevelsback+p.numinputs)), 
    rand(1:(p.numlevelsback+p.numinputs)) ]
end

# Random circuit as a list of instruction integers
# Assumes gates are arity 2
function rand_lcircuit( n_instructions::Int64, numregisters::Int64, numinputs::Int64, funcs::Vector{Func} )
  p = Parameters( numinputs, 1, n_instructions, numregisters )
  LinCircuit( [ rand_ivect( numregisters, numinputs, funcs ) for _ = 1:n_instructions ], p )
end

# Random circuit as a list of instruction integers
# Assumes gates are arity 2
# p.numinteriors is equivalent to numinstructions
# p.numlevelsback is equivalent to numregisters
function rand_lcircuit( p::Parameters, funcs::Vector{Func}=lin_funcs(p.numinputs) )
  LinCircuit( [ rand_ivect( p.numlevelsback, p.numinputs, funcs ) for _ = 1:p.numinteriors], p )
end

function execute_random_circuits( ncircuits::Int64, prog_length::Int64, numregisters::Int64, numinputs::Int64, funcs::Vector{Func} )
  n_instructions = num_instructions( numregisters, numinputs, funcs ) 
  for i = 1:ncircuits
    circuit_ints = rand(0:(n_instructions-1),prog_length)
    execute_lcircuit( circuit_ints, numregisters, numinputs, funcs )  
  end
end

# Mutates the instruction instruction_vect, returns the instruction vector of the mutated instruction
# mut_loc is the mutate location.
# Returns the mutated instruction as a vector of MyInts 
# The mutated instruction is guaranteed to be different than instruction_vect.
function mutate_instruction( instruction_vect::Vector{MyInt}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}, mut_loc::Int64; 
      nodearity::Int64=2 )  
  lfuncs = length(funcs)
  result_vect = deepcopy(instruction_vect)
  num_mutate_locations = ((lfuncs>1) ? (lfuncs-1) : 0) + (numregisters-1) + nodearity*(numregisters+numinputs-1)
  @assert 1 <= mut_loc && mut_loc <= num_mutate_locations
  if mut_loc <= lfuncs-1
    result_vect[1] = (mut_loc<instruction_vect[1]) ? mut_loc : (mut_loc+1) 
    #println("result_vect[1]: ",result_vect[1])
  elseif mut_loc <= (lfuncs-1) + (numregisters-1)
    loc = mut_loc - (lfuncs-1)
    result_vect[2] = (loc<instruction_vect[2]) ? loc : (loc+1)
    #println("result_vect[2]: ",result_vect[2])
  else
    for i = 1:nodearity
      if mut_loc <= (lfuncs-1) + (numregisters-1) + (numregisters+numinputs-1)*i
        loc = mut_loc - (lfuncs-1) - (numregisters-1) - (numregisters+numinputs-1)*(i-1)
        result_vect[2+i] = (loc<instruction_vect[2+i]) ? loc : (loc+1)
        #println("result_vect[",2+i,"]: ",result_vect[2+i])
        break
      end
    end
  end
  result_vect
end

function mutate_instruction( instruction_vect::Vector{MyInt}, p::Parameters, funcs::Vector{Func}, mut_loc::Int64; 
      nodearity::Int64=2 )  
  mutate_instruction( instruction_vect, p.numlevelsback, p.numinputs, funcs, mut_loc, nodearity=p.nodearity )
end

# Mutates the instruction instruction_vect at a random location, returns the instruction vector of the mutated instruction 
# Returns the mutated instruction as a vector of MyInts
function mutate_instruction( instruction_vect::Vector{MyInt}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )  
  lfuncs = length(funcs)
  num_mutate_locations = ((lfuncs>1) ? (lfuncs-1) : 0) + (numregisters-1) + nodearity*(numregisters+numinputs-1)
  #println("num_mutate_locations: ",num_mutate_locations)
  mut_loc = rand(1:num_mutate_locations)
  #println("num_mutate_locations: ",num_mutate_locations,"  mut_loc: ",mut_loc)
  mutate_instruction( instruction_vect, numregisters, numinputs, funcs, mut_loc, nodearity=nodearity )
end  

function mutate_instruction( instruction_vect::Vector{MyInt}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )  
  mutate_instruction( instruction_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

# Returns the vector of all possilbe mutated instructions
function mutate_instruction_all( instruction_vect::Vector{MyInt}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; 
      nodearity::Int64=2 )  
  lenfuncs = length(funcs)
  num_mutate_locations = ((lenfuncs>1) ? (lenfuncs-1) : 0) + (numregisters-1) + nodearity*(numregisters+numinputs-1)
  instructions = Vector{MyInt}[]
  for mut_loc = 1:num_mutate_locations
    new_instruction = mutate_instruction( instruction_vect, numregisters, numinputs, funcs, mut_loc, nodearity=nodearity )
    push!(instructions,new_instruction)
  end
  instructions
end

# Returns the vector of all possilbe mutated instructions
function mutate_instruction_all( instruction_vect::Vector{MyInt}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )  
  mutate_instruction_all( instruction_vect, p.numlevelsback, p.numinputs, funcs )
end

# Mutates a random instruction of circuit_vect.
function mutate_circuit!( circuit_vect::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )  
  numinstructions = length(circuit_vect)
  ind = rand(1:numinstructions)
  #println("ind: ",ind)
  circuit_vect[ind] = mutate_instruction( circuit_vect[ind], numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity=nodearity )
  circuit_vect
end 

function mutate_circuit!( circuit_vect::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )  
  mutate_circuit!( circuit_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

# Mutates a random instruction of the circuit_vect of circuit.
function mutate_circuit!( circuit::LinCircuit, funcs::Vector{Func}; nodearity::Int64=2 )
  #LinCircuit( params, mutate_circuit!( circuit.circuit_vects, circuit.params, funcs ) )
  LinCircuit( mutate_circuit!( circuit.circuit_vects, circuit.params, funcs ), circuit.params )
end

# Do all possible mutations to circuit_vect
function mutate_circuit_all( circuit_vect::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 ) 
  mutate_circuit_all( circuit_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

# Do all possible mutations to the circuit_vect of circuit
function mutate_circuit_all( circuit::LinCircuit, funcs::Vector{Func}=circuit.params.numinputs; nodearity::Int64=2 ) 
  mutate_circuit_all( circuit.circuit_vects, circuit.params, funcs )
end

# Do all possible neutral mutations to the circuit_vect of circuit
function mutate_circuit_neutral_all( circuit::LinCircuit, funcs::Vector{Func}=circuit.params.numinputs; nodearity::Int64=2 ) 
  ov = output_values( circuit )[1]
  filter( c->output_values(c)[1]==ov[1], mutate_circuit_all( circuit.circuit_vects, circuit.params, funcs ))
end

# Do all possible neutral mutations to the circuit_vect of circuit, return the circuit_ints
function mutate_circuit_neutral_all_ints( circuit::LinCircuit, funcs::Vector{Func}=circuit.params.numinputs; nodearity::Int64=2 ) 
  map( c->circuit_to_circuit_int(c,funcs), mutate_circuit_neutral_all( circuit, funcs ) )
end

# The set components is a set of LinCircuits.
# If this function doesn't stack overflow, it returns the neutral component of the neutral set of circuit.
# The LinCircuit circuit is added to components along with all mutations of circuit.
# If this did not increase the size of components, quit.
# Otherwise, recursively call this function on all of these mutations.
function neutral_component_ints!( components::Set, circuit::LinCircuit  )
  size0 = length( components )
  #println("size0: ",size0)
  mcn = mutate_circuit_neutral_all( circuit, funcs )
  mcn_ints = map( c->circuit_to_circuit_int( c, funcs ), mcn )
  union!( components, Set( mcn_ints ) )
  new_size = length( components )
  if size0 == new_size
    return components
  else
    for c in mcn
      neutral_component_ints!( components, c )
    end
  end
end

function mutate_circuit_all( circuit_vect::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 ) 
  numinstructions = length(circuit_vect)
  p = Parameters( numinputs, 1, numinstructions, numregisters )  # Assumes 1 output
  lenfuncs = length(funcs)
  #result = Vector{Vector{MyInt}}[]
  result = LinCircuit[]
  for ind = 1:length(circuit_vect)
    mutated_instructions=mutate_instruction_all(circuit_vect[ind],numregisters,numinputs,funcs,nodearity=nodearity) 
    for mi in mutated_instructions
      mutated_circuit = deepcopy(circuit_vect)
      mutated_circuit[ind] = mi
      push!(result,LinCircuit(mutated_circuit,p))
    end
  end
  result
end

# Mutates ciruit in all possible ways.
# The result depends on the keyword arguments.
# If !output_outputs && output_circuits returns list of the chromsomes produced by mutation
# If output_outputs && !output_circuits returns list of the outputs of chromsomes produced by mutation
# If output_outputs && output_circuits returns both the list of outputs and the list of chromosomes as a pair
# If robustness_only==true returns the pair: (avg_robustness, evolvability)
# Deterministic if all functions in default_funcs() have the same arity                 
function mutate_all( circuit::LinCircuit, funcs::Vector{Func}=default_funcs(circuit.params.numinputs); 
      robustness_only::Bool=false, output_outputs::Bool=true, output_circuits::Bool=false, nodearity::Int64=2 )
  #println("mutate_all: circuit: ",circuit)
  phenotype = output_values(circuit)   
  mca = mutate_circuit_all( circuit, funcs, nodearity=nodearity ) 
  if output_outputs || robustness_only
    outputs = map( x->output_values(x), mca )
  end
  if robustness_only
    output_outputs = output_circuits = false
    robustness_sum = 0.0; robustness_count = 0
    robustness_avg = length(filter(x->x==phenotype,outputs))/length(outputs)
    evolvability = length(unique(outputs))/length(outputs)
    result = (robustness_avg,evolvability)
  elseif output_outputs && !output_circuits
    result = outputs
    return result
  elseif !output_outputs && output_circuits
    result = mca
    return result
  elseif output_outputs && output_circuits
    result = (outputs,mca)
  end
end

function print_circuit( f::IO, circuit::LinCircuit, funcs::Vector{Func} )
  print(f,"(")
  i = 1
  for lc in circuit.circuit_vects
    print(f,"R$(lc[1])=$(funcs[lc[2]].name)(R$(lc[3]),R$(lc[4]))")
    if i < length(circuit.circuit_vects) print(f,",") end
    i += 1
  end
  println(f,")")
end
 
function print_circuit( circuit::LinCircuit, funcs::Vector{Func} )
  print_circuit( Base.stdout, circuit, funcs )
end

function instruction_distance( instruction_vect1::Vector{MyInt}, instruction_vect2::Vector{MyInt} )
  @assert length(instruction_vect1) == length(instruction_vect2)
  count_diffs = 0
  for i in 1:length(instruction_vect1)
    count_diffs += instruction_vect1[i] == instruction_vect2[i] ? 0 : 1
  end
  count_diffs
end

function circuit_distance( circuit_vect1::Vector{Vector{MyInt}},circuit_vect2::Vector{Vector{MyInt}} )
  @assert length(circuit_vect1) == length(circuit_vect2)
  sum_distances = 0
  for i in 1:length(circuit_vect1)
    sum_distances += instruction_distance( circuit_vect1[i], circuit_vect2[i] )
  end
  sum_distances
end  

function circuit_distance( circuit1::LinCircuit, circuit2::LinCircuit )
  circuit_distance( circuit1.circuit_vects, circuit2.circuit_vects )
end

# Returns the number of genotypes mapping to each phenotype for the parameters/funcs setting.
# This is for LinCircuits.
function count_genotypes_lc_mt( p::Parameters, funcs::Vector{Func}=lin_funcs(p.numinputs) )
  genotype_counts = [ Threads.Atomic{Int64}(0) for i= 1:2^2^p.numinputs ]
  num_circuits = count_circuits_lc( p )
  Threads.@threads for i = 1:num_circuits
    ov = output_values( circuit_int_to_circuit(i,p,funcs))[1]
    Threads.atomic_add!( genotype_counts[ov+1], 1 )
  end
  genotype_counts = map( ph->ph[], genotype_counts )
  count_genotypes_table( p, funcs, genotype_counts )   # defined in Chromosome.jl
end

count_genotype_counts_lc_mt = count_genotypes_lc_mt

# Returns the number of genotypes mapping to each phenotype for the parameters/funcs setting.
# This is for LinCircuits.
function count_genotypes_lc( p::Parameters, funcs::Vector{Func}=lin_funcs(p.numinputs) )::Vector{Int64}
  num_circuits = count_circuits_lc( p )
  genotype_counts = zeros( Int64, 2^2^p.numinputs )
  for i = 1:num_circuits
    ov = output_values( circuit_int_to_circuit(i,p,funcs))[1]
    genotype_counts[ov+1] += 1
  end
  genotype_counts
  count_genotypes_table( p, funcs, genotype_counts )
end

# Number of circuits with the given parameter setting and number of funcs
# Example:  sp = Parameters( 2, 1, 4, 2)
# count_circuits_lc( sp, nfuncs=4 ):  268435456 == 2^28 which agrees with Hu, Payne, Banzhaf 2012
function count_circuits_lc( p::Parameters; nfuncs::Int64=length(lin_funcs(p.numinputs) ) )::Int64
  @assert p.numoutputs == 1   # Not tested for more than 1 output, but probably works in this case.
  numcircuits_per_instruction = nfuncs*p.numlevelsback*(p.numlevelsback+p.numinputs)^2
  # p.numinteriors is number of instructions
  numcircuits_per_instruction^p.numinteriors^p.numoutputs
end

# Another version of count_circuits_lc()
function count_circuits_lc( p::Parameters, funcs::Vector{Func} )::Int64
  count_circuits_lc( p, nfuncs=length(funcs) )
end

# 3/1/24:  For sp = Parameters(2,1,2,2), returned 16384 circuit_vects, and count_circuits_lc(sp) == 16384
function enumerate_circuits_lc( p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs))::Vector{LinCircuit}
  ecl = enumerate_circuits_lc( p, p.numinteriors, funcs )
  map( ec->LinCircuit( ec, p ), ecl )
end

# Recursive helper function
# Returns only a circuit_vect which will be included in a LinCircuit by the top-level caller
function enumerate_circuits_lc( p::Parameters, numinstructions::Int64, funcs::Vector{Func}; maxarity::Int64=2 )
  #println("ec_lc: numinstructions: ",numinstructions)
  if numinstructions == 0
    result =  [[ MyInt[] ]]
    return result
  end
  prev_result = enumerate_circuits_lc( p, numinstructions-1, funcs )
  #println("prev_result: ",prev_result)
  new_result = Vector{Vector{MyInt}}[]
  for prev_circ in prev_result
    if length(prev_circ[1])==0
      pop!(prev_circ)
    end
    for i = 1:length(funcs)
      for j = 1:p.numlevelsback
        for h = 1:p.numlevelsback+p.numinputs
          for k = 1:p.numlevelsback+p.numinputs
            pc = push!(deepcopy(prev_circ),MyInt[i,j,h,k] )
            push!( new_result, pc )
          end
        end
      end
    end
  end
  new_result
end

