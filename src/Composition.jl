#=
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

function comp_experiment( numinputs::Int64, numcompositions::Int64, needs_list::Vector{Goal} )
  funcs = default_funcs(numinputs)
  needs = [ Need( nd, 0, Int64[] ) for nd in needs_list ]
  sort!(needs,lt=(x,y)->x.goal[1]<y.goal[1])  # assumes single component goals
  all_circuits = CompositionalCircuit[]
  for ff in funcs
    push!(all_circuits,CompGate([1,2],ff,1))
  end
  numprimitives = length(funcs)
  for i = 1:numcompositions
    ncircuits = 2
    inds = [rand(1:length(all_circuits) for _=1:ncircuits]
    inlist = vcat(collect(1:(numinputs+1)),
    new_circ = CompCircuit([rand(inlist) for _=1:numinputs],[all_circuits[ind] for in in inds],1)
  end
end

function check_need!( needs::Vector{Need}, all_circuits::Vector{CompositionalCircuit}, circ_index::Int64 )
  circ = all_circuits[circ_index]
  @assert circ.numoutputs == 1
  numinputs = length(circ.inputs)
  perms = permutations(circ.inputs)
  for p in perms
    pcirc = deepcopy(all_circuits[circ_index])
    pcirc.inputs = p
    result = execute(pcirc,construct_context(numinputs))
    need_result = Need([result],0,Int64[])
    index = searchsortedfirst(needs,need_result,lt=(x,y)->x.goal[1]<y.goal[1])
    println("perm: ",p,"  index: ",index)
    if index <= length(needs) && needs[index].goal[1] == result && needs[index].circuit_index==0  # Found an unmet need
      # Should check if number gates of circ < all_circuits[needs[index][2]] and replace if true
      println("need index: ",index)
      needs[index].circuit_index = circ_index
      needs[index].inputs = p
      break
    end
  end
end
    
function execute( circuit::CompGate, input_values::Vector{MyInt} )
  @assert circuit.numoutputs==1
  args = input_values[circuit.inputs]
  circuit.func.func(args...)
end

function execute( circuit::CompCircuit, input_values::Vector{MyInt} )
  @assert circuit.numoutputs==1
  println("circuit.circuits: ",circuit.circuits)
  extended_inputs = vcat(input_values,zeros(MyInt,length(circuit.circuits)))
  println("extended_inputs: ",extended_inputs)
  for i = 1:length(circuit.circuits)
    extended_inputs[i+length(input_values)] = execute(circuit.circuits[i],extended_inputs)
    println("extended_inputs: ",extended_inputs)
  end
  extended_inputs[end]
end

function compose_circuits( c1::CompositionalCircuit, c2::CompositionalCircuit )


function comp_circuit( ch::Chromosome )
  CompCircuit( collect(1:ch.params.numinputs), [CompGate(int.inputs,int.func,1) for int in ch.interiors], 1)
end

function print_gate_list( f::IO, node_vect::Vector{InteriorNode}, numinputs::Int64 )
  print(f,"[")
  len = length(node_vect)
  for i in 1:len
    print(f," CompGate(",node_vect[i].inputs,",",node_vect[i].func.name,",1)")
  end
  print(f,"]")
end

function print_comp_circuit( f::IO, ch::Chromosome )
  print(f,"CompCircuit(",collect(1:ch.params.numinputs),",")
  gate_list(f,ch.interiors,ch.params.numinputs)
  println(f,",1)")
end


