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
export LinCircuit, output_values
export execute_lcircuit, instruction_vect_to_instruction_int, instruction_int_to_instruction_vect, lincomplexity
export vect_to_int, int_to_vect, rand_ivect, execute_random_circuit, number_active, degeneracy, recover_phenotype
export num_instructions # This is the total number of possible instructions for a parameter setting
export rand_lcircuit, mutate_circuit!, mutate_instruction, mutate_circuit_all, mutate_all, print_circuit
export instruction_vects_to_instruction_ints, instruction_ints_to_instruction_vects
export instruction_ints_to_circuit_int, circuit_int_to_instruction_ints
export circuit_int_to_circuit, circuit_to_circuit_int, enumerate_circuits_lc, count_circuits_lc
#export circuit_int
OutputType = Int64

#=
mutable struct LinCircuit
    # Here, numinteriors represents number of gates, numlevelsback represents number of computational registers
    # the First numoutputs registers are the output registers
    circuit_vects::Vector{Vector{MyInt}}
    params::Parameters   
end
=#

function output_values( c::LinCircuit, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = lin_funcs( c.params.numinputs )
  end
  #println("funcs: ",funcs)
  R = execute_lcircuit( c, funcs )
  R[1:c.params.numoutputs]
end  

function output_values( cv::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func} )
  output_values( LinCircuit( cv, p ), funcs )
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

function lincomplexity( circuit_vects::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  p = Parameters( numinputs, 1, length(circuit_vects), numregisters )
  lincomplexity( LinCircuit(circuit_vects,p), funcs, nodearity=nodearity )
end

# Executes the circuit, returns the array X of instruction outputs
# The outputs are  R[1:numoutputs]
function lincomplexity( circuit::LinCircuit, funcs::Vector{Func}; nodearity::Int64=2 )
  p = circuit.params
  R = fill(MyInt(0), p.numlevelsback+p.numinputs )
  X = fill(MyInt(0), length(circuit.circuit_vects) )  
  R[p.numlevelsback+1:end] = construct_context(p.numinputs)
  i = 1
  for lc in circuit.circuit_vects
    #println("lc: ",lc,"  func: ",funcs[lc[1]],"  R: ",R)
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

#  Assumes 1 output registers 
#  Not finished and not debugged
#=
function number_active( circ::LinCircuit )
  p = circ.params
  @assert p.numoutputs = 1
  numoutputs = p.numoutputs
  output_instructions = Int64[]
  output_registers = Set(collect(MyInt(1):MyInt(numoutputs)))
  output_registers_found = Set(MyInt[])
  input_registers_found = Set(MyInt[])
  for i = p.numinteriors:-1:1
    println("i: ",i,"  cv: ",circ.circuit_vects[i][2],"  oi: ",output_instructions,"  or: ",output_registers,"  orf: ",output_registers_found)
    if circ.circuit_vects[i][2] in output_registers && !(circ.circuit_vects[i][2] in output_registers_found)
      push!( output_instructions, i )
      push!( output_registers_found,circ.circuit_vects[i][2] )
    end
    println("i: ",i,"  cv: ",circ.circuit_vects[i][2],"  oi: ",output_instructions,"  or: ",output_registers,"  orf: ",output_registers_found)
    output_registers = output_registers_found
    output_registers_found = Set(MyInt[])
  end
  ( output_instructions, output_registers_found )
end
=#
      
function degeneracy( circuit::LinCircuit )  # Dummy function for now 10/16/21
  return 0
end

function recover_phenotype( circuit::LinCircuit, max_steps_recover, maxtrials_recover, maxtries_recover )  # Dummy function for now 10/16/21
  return (0,0)
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
    funcs = lin_funcs( circ.params.numinputs )
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
    funcs = lin_funcs( p.numinputs )
  end
  instruction_ints_to_instruction_vects( c_ints, p.numlevelsback, p.numinputs, funcs )
end

#=
function circuit_vect_to_circuit_ints( c::LinCircuit, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = lin_funcs( c.params.numinputs )
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

function circuit_int_to_instruction_ints( c_int::Int128, p::Parameters, funcs::Vector{Func} ) 
  numinstructions = p.numinteriors
  multiplier = num_instructions( p, funcs )
  #println("multiplier: ",multiplier)
  int_to_vect( c_int, fill( multiplier, numinstructions ))
end

function circuit_to_circuit_int( circ::LinCircuit, funcs::Vector{Func} )
  instruction_ints_to_circuit_int( instruction_vects_to_instruction_ints( circ, funcs ), circ.params, funcs )
end

function circuit_int_to_circuit( c_int::Int128, p::Parameters, funcs::Vector{Func} )
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
  mut_loc = rand(1:num_mutate_locations)
  #println("num_mutate_locations: ",num_mutate_locations,"  mut_loc: ",mut_loc)
  mutate_instruction( instruction_vect, numregisters, numinputs, funcs, mut_loc, nodearity=nodearity )
end  

function mutate_instruction( instruction_vect::Vector{MyInt}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )  
  mutate_instruction( instruction_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

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

function mutate_circuit!( circuit::LinCircuit, funcs::Vector{Func}; nodearity::Int64=2 )
  #LinCircuit( params, mutate_circuit!( circuit.circuit_vects, circuit.params, funcs ) )
  LinCircuit( mutate_circuit!( circuit.circuit_vects, circuit.params, funcs ), circuit.params )
end

function mutate_circuit_all( circuit_vect::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 ) 
  mutate_circuit_all( circuit_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

function mutate_circuit_all( circuit::LinCircuit, funcs::Vector{Func}=circuit.params.numinputs; nodearity::Int64=2 ) 
  mutate_circuit_all( circuit.circuit_vects, circuit.params, funcs )
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

function count_circuits_lc( p::Parameters; nfuncs::Int64=length(default_funcs(p.numinputs) ) )
  @assert p.numoutputs == 1   # Not tested for more than 1 output, but probably works in this case.
  numcircuits_per_instruction = nfuncs*p.numlevelsback*(p.numlevelsback+p.numinputs)^2
  numcircuits_per_instruction^p.numinteriors^p.numoutputs
end

function enumerate_circuits_lc( p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs))
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

#= Moved into fnc.jl 12/28/21
function find_neutral_comps( lc_list::Vector{LinCircuit}, phenotype::MyInt, funcs::Vector{Func} )
  S = Dict{Int64,Set{Int128}}()
  p = lc_list[1].params
  new_key = 1
  for lc in lc_list
    #println("lc: ",lc)
    lci = circuit_to_circuit_int(lc,funcs)
    # mutate_all() keyword arguments don't work
    cv_list = filter( x->output_values(x,p,funcs)[1]==phenotype, mutate_circuit_all(lc.circuit_vects, p, funcs))
    #mlist = map(ci->LinCircuit(ci,p) for ci in cv_list )
    mlist = [ LinCircuit(ci,p) for ci in cv_list ]
    ihlist = map(h->circuit_to_circuit_int(h,funcs),mlist)
    push!(ihlist,circuit_to_circuit_int(lc,funcs))
    push!(ihlist,lci)
    ihset = Set(ihlist)
    if length(ihset) > 0
      for ky in keys(S)
        #println("ky: ",ky,"  S[ky]: ",S[ky])
        if length( intersect( ihset, S[ky] ) ) > 0
          union!( ihset, S[ky] )
          delete!( S, ky )
        end
      end
      S[ new_key ] = ihset
      if new_key % 100 == 0
        print("lci: ",lci,"  length(ihset): ",length(ihset),"   ")
        println("length(S[",new_key,"]) = ",length(S[new_key]))
      end
      new_key += 1
    end
  end
  for ky0 in keys(S)
    for ky1 in keys(S)
      if length(intersect(S[ky0],S[ky1]))>0 && ky0 != ky1
        println("the intersection of set S[",ky0,"] and set S[",ky1,"] is nonempty")
      end
    end
  end
  S
end  
=#
