using Permutations
using Test
using DataStructures
using Main.CGP

# The next 3 constructor functions augment the type definitions above
function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64}, numoutputs::Int64 )
  CGP.Cc( circuits, inputs, numoutputs )
end

function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64} )
  CGP.Cc( circuits, inputs, 1 )
end

function CGP.Cg( func::Func, inputs::Vector{Int64} )
  CGP.Cg( func, inputs, 1 )
end

# Runs a series of circuit compositions with a common needs list and a common all_circuits list. 
# The needs list is a list of Needs where each Need (see above) consists of a goal, a corresponding
#    circuit that maps to the goal, and the input permutation of the circuit that maps to the goal.
# Initially, needs have no corresponding circuit.
# The all_circuits list is a list of the circuits that satisfy needs and that can be combined to 
#   create new circuits.  
# It is initialized to gate circuits that implement the functions in funcs.
# needs_list is a list of (goal,numinputs) pairs where goal is the desired phenotype and numinputs is the 
#   number of inputs of the desired circuit.
# Example:  needs_list = [([0x000a],2),([0x0003],2),([0x0008],2),([0x000e],2)]
function comp_experiment( numinputs::Int64, numcompositions::Int64, needs_list::Vector{Tuple{Goal,Int64}}, 
    numcircuits::Int64, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  needs = [ CGP.Need( nd[1], 0, zeros(Int64,nd[2]) ) for nd in needs_list ]
  sort!(needs,lt=(x,y)->x.goal[1]<y.goal[1])  # assumes single component goals
  all_circuits = CGP.CC[]
  for ff in funcs
    push!(all_circuits,CGP.Cg(ff,[1,2],1))
  end
  println("all_circuits: ",all_circuits)
  for i = 1:numcompositions
    new_circ = random_comp_circuit( all_circuits, numinputs, numcircuits )
    println("ce new_circ: ",new_circ,"  result: ",execute!(new_circ,construct_context(numinputs)))
    index = check_need!(needs,all_circuits,new_circ) # also does push!(all_circuits,new_circ)
    index = 0
    if index > 0
      #push!(all_circuits,new_circ)
      println("needs updated need[",index,"]: ",needs[index])
      new_circ.inputs = needs[index].inputs
      @test execute(new_circ,construct_context(length(needs[index].inputs))) == needs[index].goal[1]
    end  
  end
  (needs, all_circuits)
end

# See function comp_experiment() for an example call with construction of all_circuits
function random_comp_circuit( all_circuits::Vector{CGP.CC}, numinputs::Int64, ncircuits::Int64,
    funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  new_circuits = CGP.CC[]
  for i = 1:ncircuits
    new_circ = deepcopy(rand(all_circuits))
    len = length(new_circ.inputs)
    #println("len: ",len,"  rc new_circ: ",new_circ)
    for j = 1:len
      new_circ.inputs[j] = rand(1:len+i-1)
    end
    push!(new_circuits,new_circ)
  end
  #println("rc new_circuits: ",new_circuits)
  CGP.CC(new_circuits,RandomPermutation(numinputs).data,1)
end

# Random CompCircuit wiith only gates as subcircuits
function random_gate_circuit( numinputs::Int64, ngates::Int64, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  all_circuits = CGP.CC[]
  for ff in funcs
    push!(all_circuits,CGP.Cg(ff,[1,2],1))
  end    
  gate_circuits = CGP.CC[]
  inlist = collect(1:numinputs)
  for j = 1:ngates
    println("j: ",j,"  numinputs+j-1: ",numinputs+j-1)
    func = funcs[rand(1:length(funcs))]
    push!(gate_circuits,CGP.Cg(func,[rand(1:(numinputs+j-1)) for _=1:func.arity],1))
  end
  new_circ = CGP.Cc(gate_circuits,RandomPermutation(numinputs).data,1)
end

#=
function check_need!( needs::Vector{CGP.Need}, all_circuits::Vector{CGP.CC}, circ::CGP.CC )
  println("circ: ",circ)
  @assert circ.numoutputs == 1
  numinputs = length(circ.inputs)
  perms = permutations(circ.inputs)
  for p in perms
    pcirc = deepcopy(circ)
    pcirc.inputs = p
    result = execute!(pcirc,construct_context(numinputs))
    need_result = CGP.Need([result],0,Int64[])
    index = searchsortedfirst(needs,need_result,lt=(x,y)->x.goal[1]<y.goal[1])
    #println("perm: ",p,"  index: ",index)
    if index <= length(needs) && needs[index].goal[1] == result && needs[index].circuit_index==0  # Found an unmet need
      # Should check if number gates of circ < all_circuits[needs[index][2]] and replace if true
      println("need index: ",index)
      push!(all_circuits,circ)
      needs[index].circuit_index = length(all_circuits)+1   # Assumes circ will be pushed onto all_circuits
      needs[index].inputs = p
      return index  # Successfully found a need
    end
  end
  return 0
end
=#

function check_need!( needs::Vector{CGP.Need}, all_circuits::Vector{CGP.CC}, circ::CGP.CC )
  #circ = all_circuits[circ_index]
  println("circ: ",circ)
  @assert circ.numoutputs == 1
  numinputs = length(circ.inputs)
  perms = permutations(circ.inputs)
  for p in perms
    pcirc = deepcopy(circ)
    pcirc.inputs = p
    result = execute!(pcirc,construct_context(numinputs))
    need_result = CGP.Need([result],0,Int64[])
    index = searchsortedfirst(needs,need_result,lt=(x,y)->x.goal[1]<y.goal[1])
    #println("perm: ",p,"  index: ",index)
    if index <= length(needs) && needs[index].goal[1] == result && needs[index].circuit_index==0  # Found an unmet need
      # Should check if number gates of circ < all_circuits[needs[index][2]] and replace if true
      println("need index: ",index)
      push!(all_circuits,circ)
      needs[index].circuit_index = length(all_circuits)+1   # Assumes circ will be pushed onto all_circuits
      needs[index].inputs = p
      return index  # Successfully found a need
    end
  end
  return 0
end

# Execute CompGate circuit and return values when circuit is applied to context (which may be input_values).
# extended_inputs contains the values computed by earlier circuits.
function execute!( circuit::CGP.Cg, input_values::Vector{MyInt}, extended_inputs::Vector{MyInt}=MyInt[]; 
     active_only::Bool=true, use_inputs::Bool=true  )
  @assert circuit.numoutputs==1
  println("gate execute circuit: ",circuit,"  extended inputs: ",extended_inputs)
  if length(extended_inputs) == 0
    extended_inputs = input_values 
  end
  println("gate extended_inputs: ",extended_inputs)
  args = use_inputs ? extended_inputs[circuit.inputs] : extended_inputs
  println("gate execute args: ",args)
  circuit.func.func(args...)
end

# Execute CompCircuit circuit and return values when circuit is applied to input_values (which is the context on the original call).
# extended_inputs contains the values computed by earlier circuits.
# If use_inputs==true, then the inputs components of the circuit is used to reorder the input_values components.
# If active_only==false, starts with the last circuit and recursively computes the values needed.
# If active_only==true, starts with the input_values and iterates through the circuits, but still uses recursion.
function execute!( circuit::CGP.CC, input_values::Vector{MyInt}, extended_inputs::Vector{MyInt}=MyInt[]; 
    use_inputs::Bool=true, active_only::Bool=true)
  @assert circuit.numoutputs==1
  println("circuit execute circuit: ",circuit)
  if length(extended_inputs) == 0
    extended_inputs = vcat( input_values, zeros(MyInt,length(circuit.circuits)) )
  end
  args = use_inputs ? extended_inputs[circuit.inputs] : extended_imputs
  extended_inputs[1:length(args)] = args
  println("args: ",args,"  extended_inputs: ",extended_inputs)
  if active_only
    for i in active_list( circuit ) 
      if i > length(args)
        ext_inputs = typeof((circuit.circuits[i-length(args)])) == CGP.Cg ? extended_inputs : MyInt[]
        extended_inputs[i] = execute!(circuit.circuits[i-length(args)],args,ext_inputs,active_only=active_only)
      end
      println("i: ",i,"  extended_inputs: ",extended_inputs)
    end
  else
    for i = 1:length(circuit.circuits)
      ext_inputs = typeof(circuit.circuits[i]) == CGP.Cg ? extended_inputs : MyInt[]
      extended_inputs[i+length(args)] = execute!(circuit.circuits[i],args,ext_inputs,active_only=active_only)
      println("i: ",i,"  extended_inputs: ",extended_inputs)
    end
  end
  println("return value: ",@sprintf("0x00%x",extended_inputs[end]))
  extended_inputs[end]
end

# computes the indices of the circuit components that need to be coputed to compute the value of the given circuit
function active_list( circuit::CGP.CC )
  input_set = Set(Int64[])
  for index = length(circuit.circuits):-1:1
    for i = 1:length(circuit.circuits[index].inputs)
      input_set = union!(input_set,circuit.circuits[index].inputs[i])
    end
  end
  input_set = union!(input_set,length(circuit.inputs)+length(circuit.circuits))  # Last circuit is always active
  input_indices = sort([i for i in input_set])
end

# Convert a chromsome to a CompCircuit
# The inputs of the CompCircuit are always [1,2....]
function chromosome_to_circuit( ch::Chromosome )
  inputs = [ inp.index for inp in ch.inputs ]
  gates = CGP.CC[ CGP.Cg( i.func, i.inputs, 1 ) for i in ch.interiors ]
  numoutputs = length(ch.outputs)
  CGP.CC( gates, inputs, numoutputs )
end

function circuit_to_chromosome( circuit::CGP.CC )
  int_list = InteriorNode[]
  params = Parameters( length(circuit.inputs), 1, length(circuit.circuits), length(circuit.circuits))
  input_list = [ InputNode( i ) for i in circuit.inputs ]
  last_circuit = length(circuit.inputs)+length(circuit.circuits)
  first_output = last_circuit-circuit.numoutputs+1
  output_range = first_output:last_circuit
  output_list = [ OutputNode(i) for i = output_range ]
  for c in circuit.circuits
    @assert typeof(c) == CGP.Cg
    int = InteriorNode( c.func, c.inputs, false, MyInt(0) )
    push!(int_list,int)
  end
  Chromosome( params, input_list, int_list, output_list, 0.0, 0.0 )
end      

# "prints" a compositional circuit to a string
function sprint_comp_circuit( cc::CGP.Cc )
  result = "CGP.CC(CGP.CC["
  for i in 1:length(cc.circuits)
    csub = cc.circuits[i]
    if typeof(csub) == CGP.Cg
      result *= sprint_gate_circuit( csub )
    elseif typeof(csub) == CGP.Cc
      result *= sprint_comp_circuit( csub )
    end
    if i < length(cc.circuits)
      result *= ", "
    end 
  end
  result *= "], " * sprint_int_list(cc.inputs) * ")"
end

function sprint_gate_circuit( cg::CGP.Cg )
  result = "CGP.Cg(" * cg.func.name * "," * sprint_int_list(cg.inputs) * ")"
end    
    
  
# Revised 8/21/21
# Converts a chromosome to a circuit and prints that circuit to IO.
# Does not print the 1 as the last component
function print_comp_circuit( f::IO, ch::Chromosome )
  print(f,"CGP.CC( ")
  print_gate_list(f,ch.interiors,ch.params.numinputs)
  # println(f,", ", sprint_int_list(collect(1:ch.params.numinputs)),", 1)")
  println(f,", ", sprint_int_list(collect(1:ch.params.numinputs)),")")
end

# Helper function for print_comp_circuit
function print_gate_list( f::IO, node_vect::Vector{InteriorNode}, numinputs::Int64 )
  print(f,"CCP.CC[")
  len = length(node_vect)
  for i in 1:(len-1)
    # print(f," CGP.Cg(",node_vect[i].func.name,",",sprint_int_list(node_vect[i].inputs),",1),")
    print(f," CGP.Cg(",node_vect[i].func.name,",",sprint_int_list(node_vect[i].inputs),"),")
  end
  # print(f," CGP.Cg(",node_vect[end].func.name,",",sprint_int_list(node_vect[end].inputs),",1)")
  print(f," CGP.Cg(",node_vect[end].func.name,",",sprint_int_list(node_vect[end].inputs),")")
  print(f,"]")
end

# String prints an Integer list with no "Integer" prefix and no extra spaces.
function sprint_int_list( lst::Vector{Union{Int64,Integer}} )
  result = "["
  for i in collect(1:(length(lst)))
    result *= @sprintf("%d",lst[i])
    if i != length(lst)
      result *= ","
    end
  end
  result *= "]"
end 
