using Permutations
using Test
using DataStructures
#=  See aliases for the loaded version
abstract type CompositionalCircuit end
mutable struct CompGate <: CompositionalCircuit
  inputs::Vector{Int64}
  func::Func
  numoutputs::Int64
end
mutable struct CompCircuit <: CompositionalCircuit
  inputs::Vector{Int64}
  circuits::Vector{CompositionalCircuit}
  numoutputs::Int64
end
mutable struct Need
  goal::Goal
  circuit_index::Int64   # If nonzero, the index of the element of all_circuits that meets the need
  inputs::Vector{Int64}
end
=#

# Runs a series of circuit compositions with a common needs list and a common all_circuits list. 
# The needs list is a list of Needs where each Need (see above) consists of a goal, a corresponding
#    circuit that maps to the goal, and the input permutation of the circuit that maps to the goal.
# Initially, needs have no corresponding circuit.
# The all_circuits list is a list of the circuits that satisfy needs and that can be combined to 
#   create new circuits.  
# It is initialized to gate circuits that implement the functions in funcs.
# needs_list is a list of (goal,numinputs) pairs where goal is the desired phenotype and numinputs is the 
#   number of inputs of the desired circuit.
function comp_experiment( numinputs::Int64, numcompositions::Int64, needs_list::Vector{Tuple{Goal,Int64}}, 
    numcircuits::Int64, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  needs = [ CGP.Need( nd[1], 0, zeros(Int64,nd[2]) ) for nd in needs_list ]
  sort!(needs,lt=(x,y)->x.goal[1]<y.goal[1])  # assumes single component goals
  all_circuits = CompositionalCircuit[]
  for ff in funcs
    push!(all_circuits,CGP.CompGate([1,2],ff,1))
  end
  println("all_circuits: ",all_circuits)
  for i = 1:numcompositions
    new_circ = random_comp_circuit( all_circuits, numinputs, numcircuits )
    println("ce new_circ: ",new_circ,"  result: ",execute(new_circ,construct_context(numinputs)))
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

function random_comp_circuit( all_circuits::Vector{CompositionalCircuit}, numinputs::Int64, ncircuits::Int64,
    funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  new_circuits = CompositionalCircuit[]
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
  CGP.CompCircuit(RandomPermutation(numinputs).data,new_circuits,1)
end

# Random CompCircuit wiith only gates as subcircuits
function random_gate_circuit( numinputs::Int64, ngates::Int64, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  all_circuits = CompositionalCircuit[]
  for ff in funcs
    push!(all_circuits,CGP.CompGate([1,2],ff,1))
  end    
  gate_circuits = CompGate[]
  inlist = collect(1:numinputs)
  for j = 1:ngates
    println("j: ",j,"  numinputs+j-1: ",numinputs+j-1)
    func = funcs[rand(1:length(funcs))]
    push!(gate_circuits,CompGate([rand(1:(numinputs+j-1)) for _=1:func.arity],func,1))
  end
  new_circ = CGP.CompCircuit(RandomPermutation(numinputs).data,gate_circuits,1)
end

function check_need!( needs::Vector{CGP.Need}, all_circuits::Vector{CompositionalCircuit}, circ::CompositionalCircuit )
  #circ = all_circuits[circ_index]
  println("circ: ",circ)
  @assert circ.numoutputs == 1
  numinputs = length(circ.inputs)
  perms = permutations(circ.inputs)
  for p in perms
    pcirc = deepcopy(circ)
    pcirc.inputs = p
    result = execute(pcirc,construct_context(numinputs))
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

function execute!( circuit::CGP.CompGate, context::Vector{MyInt}, extended_inputs::Vector{MyInt}=MyInt[]; 
    use_inputs::Bool=true  )
  @assert circuit.numoutputs==1
  println("execute circuit: ",circuit)
  if length(extended_inputs) == 0
    extended_inputs = context 
  end
  args = use_inputs ? extended_inputs[circuit.inputs] : extended_inputs
  println("execute args: ",args)
  circuit.func.func(args...)
end

# Execute circuit and return values when circuit is applied to input_values.
# input_values may be the context.
# This version starts with the inputs and iteratively computes the values of the circuits.
function execute!( circuit::CGP.CompCircuit, context::Vector{MyInt}, extended_inputs::Vector{MyInt}=MyInt[]; 
    use_inputs::Bool=true, active_only::Bool=true)
  @assert circuit.numoutputs==1
  println("execute circuit: ",circuit)
  if length(extended_inputs) == 0
    extended_inputs = vcat( context, zeros(MyInt,length(circuit.circuits)) )
  end
  args = use_inputs ? extended_inputs[circuit.inputs] : extended_imputs
  println("extended_inputs: ",extended_inputs)
  if active_only
    for i in active_list( circuit ) 
      if i > length(args)
        extended_inputs[i] = execute!(circuit.circuits[i-length(args)],extended_inputs,active_only=active_only)
      end
      println("i: ",i,"  extended_inputs: ",extended_inputs)
    end
  else
    for i = 1:length(circuit.circuits)
      extended_inputs[i+length(args)] = execute!(circuit.circuits[i],extended_inputs,active_only=active_only:)
      println("i: ",i,"  extended_inputs: ",extended_inputs)
    end
  end
  println("execute: ",@sprintf("0x00%x",extended_inputs[end]))
  extended_inputs[end]
end

function active_list( circuit::CGP.CompCircuit )
  input_set = Set(Int64[])
  for index = length(circuit.circuits):-1:1
    for i = 1:length(circuit.circuits[index].inputs)
      input_set = union!(input_set,circuit.circuits[index].inputs[i])
    end
  end
  input_set = union!(input_set,length(circuit.inputs)+length(circuit.circuits))  # Last circuit is always active
  input_indices = sort([i for i in input_set])
end

function print_comp_circuit( f::IO, ch::Chromosome )
  print(f,"CGP.CompCircuit(",collect(1:ch.params.numinputs),",")
  print_gate_list(f,ch.interiors,ch.params.numinputs)
  println(f,",1)")
end

# Helper function for print_comp_circuit
function print_gate_list( f::IO, node_vect::Vector{InteriorNode}, numinputs::Int64 )
  print(f,"[")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,"CompGate(",node_vect[i].inputs,",",node_vect[i].func.name,",1),")
  end
  print(f,"CompGate(",node_vect[len].inputs,",",node_vect[len].func.name,"1)")
  print(f,"]")
end

# Convert a chromsome to a CompCircuit
# The inputs of the CompCircuit are always [1,2....]
function comp_circuit( ch::Chromosome )
  inputs = collect(1:ch.params.numinputs)
  gates = [ CompGate( i.inputs, i.func, 1 ) for i in ch.interiors ]
  numoutputs = length(ch.outputs)
  CompCircuit( inputs, gates, numoutputs )
end

# Find those inputs from c.inputs which are actually used in the execution of c
function find_inputs( c::CompositionalCircuit )
  input_dict = Dict{Int64,Int64}()
  for cc in c.circuits
    for i in cc.inputs
      if i <= length( c.inputs )
        input_dict[i] = 0
      end
    end
  end
  kys = [ k for k in  keys(input_dict)]
  println("inputs: ",c.inputs,"  keys: ",kys)
  return c.inputs[kys]
end

# Find active subcircuits
function find_active( c::CompositionalCircuit )
  active_indices = Accumulator{Int64,Int64}()
  find_active_helper!( active_indices, [length(c.circuits)], c.circuits, 8 )
  active_indices
end
  
function find_active_helper!( active_indices::Accumulator{Int64,Int64}, cindices::Vector{Int64}, 
    clist::Vector{CompositionalCircuit}, count::Int64 )
  if count <= 0 return end
  println("cindices: ",cindices)
  for i in cindices
    println("i: ",i,"  clist[i].inputs: ",clist[i].inputs)
    if i > length(c.inputs)
      inc!( active_indices, i )
      find_active_helper!( active_indices, clist[i].inputs, clist, count-1 )
    end
  end
end
  
