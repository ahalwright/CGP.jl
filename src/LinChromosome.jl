# Linear GP circuit chromosome functions
# A circuit is a sequence of instructions on the register array R 
#   which has length numregisters + numinputs
# R is normally local to the function where is is used.
# The first numregisters components or R are calculation registers.
# The remaining numinputs components are read-only input registers which are 
#   initialized with the components of the current context.
# Even though there is a nodearity keyword value, this file assumes nodearity==2
# An instruction is specified either by a 4-tuple of MyInts or by an integer of OutputType.
# The elements of the 4-tuple are:
#  1.  The index of the logical operation (such as AND, OR) in the funcs array.
#  2.  The index of the element of R where the output is stored
#  3.  The index of the element of R which is the first operand of the logical operation
#  4.  The index of the element of R which is the second operand of the logical operation
# The functions vect_to_int() and int_to_vect() convert between these representations.
# See notes/3_24.txt for outline
# See test/testLinCircuit.jl for tests.
export LinCircuit, output_values
export execute_lcircuit, numinstructions, vect_to_int, int_to_vect, rand_ivect, execute_random_circuit
export rand_lcircuit, mutate_circuit!, mutate_instruction, mutate_all, print_circuit
OutputType = Int64

#=
mutable struct LinCircuit
    # Here, numinteriors represents number of gates, numlevelsback represents number of computational registers
    # the First numoutputs registers are the output registers
    params::Parameters   
    circuit_vects::Vector{Vector{MyInt}}
end
=#

function output_values( c::LinCircuit )
  funcs = lin_funcs( c.params.numinputs )
  #println("funcs: ",funcs)
  R = execute_lcircuit( c.circuit_vects, c.params.numlevelsback, c.params.numinputs, funcs, nodearity=c.params.nodearity )
  R[1:c.params.numoutputs]
end  

function output_values( c::LinCircuit, funcs::Vector{Func} )
  R = execute_lcircuit( c.circuit_vects, c.params.numlevelsback, c.params.numinputs, funcs, nodearity=c.params.nodearity )
  R[1:c.params.numoutputs]
end  

# Not used.  Assumes that R and funcs are in the execution environment
#function vect_to_funct( lc::Vector{MyInt}, funcs::Vector{Func} )
function vect_to_funct( lc::Vector{MyInt} )
  (lc)->R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]]) 
end

# Executes the circuit, returns the register array R.
# numregisters is the number of calculation registers, not the total number of registers
# The outputs are  R[1:numoutputs]
function execute_lcircuit( circuit_ints::Vector{OutputType}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  lcv = map( ci->int_to_vect( ci, numregisters, numinputs, funcs ), circuit_ints )
  execute_lcircuit( lcv, numregisters, numinputs, funcs )
end

#function execute_lcircuit( circuit_vects::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func} )
function execute_lcircuit( circuit::LinCircuit, p::Parameters, funcs::Vector{Func} )
  R = fill(MyInt(0), p.numlevelsback + p.numinputs )
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
    #println("lc: ",lc,"  R: ",R)
    R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
  end
  R
end

# The number of possible instructions with these settings.
function numinstructions( numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; 
     nodearity::Int64=2 )
  length(funcs)*numregisters*(numinputs+numregisters)^nodearity
end

function numinstructions( p::Parameters, funcs::Vector{Func} )
  numinstructions( p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

# Converts an instruction specified by a vector to the instruction specified by an integer.
function vect_to_int( inst_vect::Vector{MyInt}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  result = OutputType(inst_vect[1]-1)
  multiplier = 1
  #println("multiplier: ",multiplier,"  result: ",result)
  multiplier = numregisters
  result *= multiplier
  result += inst_vect[2]-1
  #println("multiplier: ",multiplier,"  result: ",result)
  #multiplier = numregisters
  multiplier = numregisters+numinputs
  for j = 1:nodearity
    result *= multiplier
    result += inst_vect[2+j]-1
    #println("multiplier: ",multiplier,"  result: ",result)
    multiplier = numregisters+numinputs
  end
  result
end

# Converts an instruction specified by an integer to the instruction specified by a vector.
function int_to_vect( inst_int::OutputType, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )  
  result = fill(MyInt(0),2+nodearity)
  i = 2 + nodearity
  for j = 1:nodearity
    multiplier = numregisters+numinputs
    result[i] = inst_int % multiplier+1
    inst_int = div(inst_int,multiplier)
    #println("i: ",i,"  multiplier: ", multiplier,"  inst_int: ",inst_int,"  result[i]: ",result[i])
    i -= 1
  end
  multiplier = numregisters
  result[i] = inst_int % multiplier+1
  inst_int = div(inst_int,multiplier)
  #println("i: ",i,"  multiplier: ", multiplier,"  inst_int: ",inst_int,"  result[i]: ",result[i])
  i -= 1
  multiplier = length(funcs)
  result[i] = inst_int % multiplier+1
  inst_int = div(inst_int,multiplier)
  #println("i: ",i,"  multiplier: ", multiplier,"  inst_int: ",inst_int,"  result[i]: ",result[i])
  result
end 
  
# Random circuit as a list of instruction vectors
# Assumes gates are arity 2
function rand_ivect( numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  MyInt[ rand(1:length(funcs)), rand(1:numregisters),  rand(1:(numregisters+numinputs)), rand(1:(numregisters+numinputs)) ]
end

# Random circuit as a list of instruction vectors
# Assumes gates are arity 2
function rand_ivect( p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )
  MyInt[ rand(1:length(funcs)), rand(1:p.numlevelsback),  rand(1:(p.numlevelsback+p.numinputs)), 
    rand(1:(p.numlevelsback+p.numinputs)) ]
end

# Random circuit as a list of instruction integers
# Assumes gates are arity 2
function rand_lcircuit( n_instructions::Int64, numregisters::Int64, numinputs::Int64, funcs::Vector{Func} )
  p = Parameters( numinputs, 1, n_instructions, numregisters )
  LinCircuit( p, [ rand_ivect( numregisters, numinputs, funcs ) for _ = 1:n_instructions ] )
end

function rand_lcircuit( p::Parameters, funcs::Vector{Func} )
  LinCircuit( p, [ rand_ivect( p.numlevelsback, p.numinputs, funcs ) for _ = 1:p.numinteriors] )
end

# Not used:  consider deleting
function lin_chromosome( lcircuits::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func} )
  R = fill(MyInt(0), numregisters+numinputs )
  R[numregisters+1:end] = construct_context(numinputs)
  #lfuncts = [  (lc)->R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]]) for ic in lcircuits ]
  println("R: ",R)
  for lc in lcircuits
    R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
    println("R: ",R)
  end
  R
end

function execute_random_circuits( ncircuits::Int64, prog_length::Int64, numregisters::Int64, numinputs::Int64, funcs::Vector{Func} )
  n_instructions = numinstructions( numregisters, numinputs, funcs ) 
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
    done = (result_vect != instruction_vect)
  result_vect
end

function mutate_instruction( instruction_vect::Vector{MyInt}, p::Parameters, funcs::Vector{Func}, mut_loc::Int64; 
      nodearity::Int64=2 )  
  mutate_instruction( instruction_vect, p.numlevelsback, p.numinputs, funcs, mut_loc, nodearity=p.nodearity )
end

# Mutates the instruction instruction_vect at a random location, returns the instruction vector of the mutated instruction 
# Returns the mutated instruction as a vector of MyInts
# There is no restriction that the mutated instruction is different from instruction_vect
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

function mutate_circuit!( instruction_vect::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )  
  numinstructions = length(instruction_vect)
  ind = rand(1:numinstructions)
  #println("ind: ",ind)
  instruction_vect[ind] = mutate_instruction( instruction_vect[ind], numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity=nodearity )
  instruction_vect
end 

function mutate_circuit!( instruction_vect::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 )  
  mutate_circuit!( instruction_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
end

function mutate_circuit!( circuit::LinCircuit, funcs::Vector{Func}; nodearity::Int64=2 )
  LinCircuit( circuit.params, mutate_circuit!( circuit.circuit_vects, circuit.params, funcs ))
end

function mutate_all( instruction_vect::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 ) 
  lfuncs = length(funcs)
  result = Vector{MyInt}[]
  num_mutate_locations = ((lfuncs>1) ? (lfuncs-1) : 0) + (numregisters-numinputs-1) + nodearity*(numregisters-1)
  for ind = 1:length(instruction_vect)
    for mut_loc = 1:num_mutate_locations
      mi = mutate_instruction(instruction_vect[ind],numregisters,numinputs,funcs,mut_loc,nodearity=nodearity)
      #println("ind: ",ind,"  mut_loc: ",mut_loc,"  mi: ",mi)
      push!(result,mi)
    end
  end
  result
end

function mutate_all( instruction_vect::Vector{Vector{MyInt}}, p::Parameters, funcs::Vector{Func}; nodearity::Int64=2 ) 
  mutate_all( instruction_vect, p.numlevelsback, p.numinputs, funcs, nodearity=p.nodearity )
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
