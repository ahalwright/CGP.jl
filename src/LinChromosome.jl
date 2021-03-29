# Linear GP circuit chromosome functions
# A circuit is a sequence of instructions on the register array R 
#   which has length numregisters + numinputs
# R is normally local to the function where is is used.
# The first numregisters components or R are calculation registers.
# The remaining numinputs components are read-only input registers which are 
#   initialized with the components of thecurrent context.
# Even though there is a nodearity keyword value, this file assumes nodearity==2
# An instruction is specified either by a 4-tuple of MyInts or by an integer of OutputType.
# The elements of the 4-tuple are:
#  1.  The index of the logical operation (such as AND, OR) in the funcs array.
#  2.  The index of the element of R where the output is stored
#  3.  The index of the element of R which is the first operand of the logical operation
#  4.  The index of the element of R which is the second operand of the logical operation
# The functions vect_to_int() and int_to_vect() convert between these representations.
# See notes/3_24.txt for outline
export execute_lcircuit, numinstructions, vect_to_int, int_to_vect, rand_ivect, execute_random_circuit
export rand_lcircuit
OutputType = Int64

# Not used.  Assumes that R and funcs are in the execution environment
#function vect_to_funct( lc::Vector{MyInt}, funcs::Vector{Func} )
function vect_to_funct( lc::Vector{MyInt} )
  (lc)->R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]]) 
end

# Executes the circuit, returns the register array R.
# The outputs are  R[1:numoutputs]
function execute_lcircuit( circuit_ints::Vector{OutputType}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  lcv = map( ci->int_to_vect( ci, numregisters, numinputs, funcs ), circuit_ints )
  execute_lcircuit( lcv, numregisters, numinputs, funcs )
end

# Executes the circuit, returns the register array R.
# The outputs are  R[1:numoutputs]
function execute_lcircuit( circuit_vects::Vector{Vector{MyInt}}, numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; nodearity::Int64=2 )
  R = fill(MyInt(0), numregisters+numinputs )
  R[numregisters+1:end] = construct_context(numinputs)
  for lc in circuit_vects
    R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
    #println("ci: ",ci,"  lc: ",lc,"  R: ",R)
  end
  R
end

# The number of possible instructions with these settings.
function numinstructions( numregisters::Int64, numinputs::Int64, funcs::Vector{Func}; 
     nodearity::Int64=2 )
  length(funcs)*numregisters*(numinputs+numregisters)^nodearity
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

# Random circuit as a list of instruction integers
# Assumes gates are arity 2
function rand_lcircuit( n_instructions::Int64, numregisters::Int64, numinputs::Int64, funcs::Vector{Func} )
  [ rand_ivect( numregisters, numinputs, funcs ) for _ = 1:n_instructions ]
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
  
