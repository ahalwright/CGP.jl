using Permutations
using Test
using DataStructures
using Main.CGP

export CC, Cg, comp_experiment, random_comp_circuit, random_gate_circuit, check_need!, execute!
export add_needs!, count_gates
export active_list, chromosome_to_circuit, circuit_to_chromosome, sprint_comp_circuit
export sprint_gate_circuit, print_comp_circuit, print_gate_list, sprint_int_list

# The next constructor functions augment the type definitions above
function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64}, numoutputs::Int64 )
  CGP.Cc( circuits, inputs, numoutputs )
end

function CGP.Cg( func::Func, inputs::Vector{Int64}, numoutputs::Int64 )
  CGP.Cg( func, inputs, numoutputs, [MyInt(0)], "" )
end

function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64} )
  CGP.Cc( circuits, inputs, 1 )
end

function CGP.Cg( func::Func, inputs::Vector{Int64} )
  CGP.Cg( func, inputs, 1 )
end

function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64}, numoutputs::Int64 )
  CGP.Cc( circuits, inputs, numoutputs, [MyInt(0)], "", 0 )
end

function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64} )
  CGP.Cc( circuits, inputs, 1, [MyInt(0)], "", 0 )
end

function CGP.CC( circuits::Vector{CC}, inputs::Vector{Int64}, cache::Goal, cache_valid::Bool )
  CGP.Cc( circuits, inputs, 1, cache, "", cache_valid )
end

function CGP.Cg( func::Func, inputs::Vector{Int64}, numoutputs::Int64 )
  CGP.Cg( func, inputs, numoutputs, [MyInt(0)], "" )
end

function CGP.Cg( func::Func, inputs::Vector{Int64} )
  CGP.Cg( func, inputs, 1, [MyInt(0)], "" )
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
# Examples:  
#    needs_list = [([0x000a],2),([0x0003],2),([0x0008],2),([0x000e],2)]
#    needs_list = [([0x000a],2),([0x0003],2),([0x0008],2),([0x000e],2),([0x0006],2),([0x0009],2)]
function comp_experiment( numinputs::Int64, numcompositions::Int64, needs_list::Vector{Tuple{Goal,Int64}},
    numcircuits::Int64, funcs::Vector{Func}=Func[]; recurse_limit::Int64=50 )
  recurse_step = Int64[0]
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  ctx = construct_context(numinputs)
  key_array = [0]
  recurse_step = [0]
  needs = Dict{Goal,Need}()
  all_circuits = Dict{String,CC}()
  i = 1
  for nd in needs_list
    needs[nd[1]] = Need(nd[1],collect(1:nd[2]),"",0.0)
    i += 1
  end
  for func in funcs
    key = func.name
    new_circ = Cg( func, collect(1:func.arity), 1, [0x0000], key )
    all_circuits[key] = new_circ
    func_goal = [CGP.execute!(all_circuits[key],ctx,recurse_step,needs,active_only=false,recurse_limit=recurse_limit)]
    need = get( needs, func_goal, nothing )
    if need != nothing
      need.circuit_key = key
      need.best_fitness = 1.0
    end
  end
  println("ce recurse_limit: ",recurse_limit)
  new_circ = Cc( CC[], collect(1:numinputs), 1, [0x0000], "", 0 )
  for i = 1:numcompositions
    rand_key = rand(keys(all_circuits))
    while rand_key == new_circ.key
      rand_key = rand(keys(all_circuits))
    end  
    println("rand_key: ",rand_key)
    new_subcirc = all_circuits[rand_key]
    new_subcirc.inputs = [ rand(1:(numinputs+i-1)) for _=1:length(new_subcirc.inputs)]
    push!( new_circ.circuits, new_subcirc )
    println("i: ",i,"  new_subcirc: ",new_subcirc)
    # also adds new_circ to all_circuits if it meets a need
    add_need!(needs,all_circuits,new_circ,key_array,recurse_step,recurse_limit=recurse_limit) 
  end
  (needs, all_circuits)
end 

function new_circuit_key( key_array::Vector{Int64} )
  key_array[1] += 1
end

# Execute CompGate circuit and return values when circuit is applied to context (which may be input_values).
# extended_inputs contains the values computed by earlier circuits.
function execute!( circuit::CGP.Cg, input_values::Vector{MyInt}, recurse_step::Vector{Int64}, needs::Dict{Array{UInt16,1},Need}, extended_inputs::Vector{MyInt}=MyInt[]; active_only::Bool=true, use_inputs::Bool=true, recurse_limit::Int64=50  )
  @assert circuit.numoutputs==1
  if length(extended_inputs) == 0
    extended_inputs = input_values 
  end
  #println("gate extended_inputs: ",extended_inputs)
  args = use_inputs ? extended_inputs[circuit.inputs] : extended_inputs
  #println("gate execute args: ",args)
  println("execute gate result: ",circuit.func.func(args...))
  circuit.func.func(args...)
end

# version of execute! without the needs parameter
function execute!( circuit::CGP.Cg, input_values::Vector{MyInt}, recurse_step::Vector{Int64}, extended_inputs::Vector{MyInt}=MyInt[]; active_only::Bool=true, use_inputs::Bool=true, recurse_limit::Int64=50  )
  needs = Dict{Goal,Need}()
  execute!( circuit, input_values, recurse_step, needs, extended_inputs, active_only=active_only, use_inputs=use_inputs, recurse_limit=recurse_limit )
end

# Execute CompCircuit circuit and return values when circuit is applied to input_values (which is the context on the original call).
# extended_inputs contains the values computed by earlier circuits.
# If use_inputs==true, then the inputs components of the circuit is used to reorder the input_values components.
# If active_only==false, starts with the last circuit and recursively computes the values needed.
# If active_only==true, starts with the input_values and iterates through the circuits, but still uses recursion.
function execute!( circuit::CGP.CC, input_values::Vector{MyInt}, recurse_step::Vector{Int64}, needs::Dict{Array{UInt16,1},Need}, extended_inputs::Vector{MyInt}=MyInt[]; use_inputs::Bool=true, active_only::Bool=true, recurse_limit::Int64=50 )
  num_gates = Int64[0]
  println("execute CC circ count gates: ",count_gates(circuit,num_gates))
  println("ex recurse_limit: ",recurse_limit)
  if circuit.cache_valid != 0
    println("ex CC returning output ",output, " as cache.")
    return circuit.output
  end
  @assert circuit.numoutputs==1
  recurse_step[1] += 1
  println("execute! recurse step incremented to ",recurse_step[1])
  if recurse_step[1] >= recurse_limit
    println("execute! exited with recurse_step ",recurse_step[1])
    return
  end
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
        println("i: ",i,"  circuit execute 187 circuit: ",circuit.circuits[i-length(args)])
        extended_inputs[i] = execute!(circuit.circuits[i-length(args)],args,recurse_step,needs,ext_inputs,active_only=active_only,recurse_limit=recurse_limit)
        println("i: ",i,"  extended_inputs: ",extended_inputs)
      end
    end
  else
    for i = 1:length(circuit.circuits)
      ext_inputs = typeof(circuit.circuits[i]) == CGP.Cg ? extended_inputs : MyInt[]
      println("i: ",i,"  circuit execute 197 circuit: ",circuit.circuits[i])
      extended_inputs[i+length(args)] = execute!(circuit.circuits[i],args,recurse_step,needs,ext_inputs,active_only=active_only,recurse_limit=recurse_limit)
      println("i: ",i,"  extended_inputs: ",extended_inputs)
    end
  end
  println("return value: ",@sprintf("0x00%x",extended_inputs[end]))
  circuit.output = [extended_inputs[end]]
  return extended_inputs[end]
end

# version of execute! without the needs parameter
function execute!( circuit::CGP.CC, input_values::Vector{MyInt}, recurse_step::Vector{Int64}, extended_inputs::Vector{MyInt}=MyInt[]; active_only::Bool=true, use_inputs::Bool=true, recurse_limit::Int64=50  )
  needs = Dict{Goal,Need}()
  execute!( circuit, input_values, recurse_step, needs, extended_inputs, active_only=active_only, use_inputs=use_inputs, recurse_limit=recurse_limit )
end

# Does not consider approximate goal matches:  best_fitnesses are set to 1.0
function add_need!( needs::Dict{Array{MyInt,1},Need}, all_circuits::Dict{String,CC}, circ::CC, key_array::Vector{Int64}, 
    recurse_step::Vector{Int64}; recurse_limit::Int64=50 ) 
  println("circ: ",circ)
  println("length(circ.circuits): ",length(circ.circuits))
  #println("length(circ.circuits[1].circuits)): ",length(circ.circuits[1].circuits))
  println("an recurse_limit: ",recurse_limit)
  recurse_step[1] += 1
  println("add_need recurse step incremented to ",recurse_step[1])
  if recurse_step[1] >= recurse_limit
    println("add_need! exited with recurse_step ",recurse_step[1])
    return
  end
  @assert circ.numoutputs == 1
  add_need = false
  numinputs = length(circ.inputs)
  output_goal = [CGP.execute!(circ,construct_context(numinputs),recurse_step,needs,active_only=false,recurse_limit=recurse_limit)]
  circ.output = output_goal
  num_gates = [0]
  new_need = Need(output_goal,circ.inputs,circ.key,1.0)
  if haskey(needs,output_goal)  # Nothing is added if false
    old_need = needs[output_goal]
    println("old_need: ",old_need)
    if old_need.best_fitness == 0.0
      #=
      if typeof(circ) == Cg
        println("pop name: ",circ.func.name)
        pop!(all_circuits,circ.func.name,nothing)
      end
      =#
      add_need = true
    elseif old_need.best_fitness < new_need.best_fitness
      add_need = true
    elseif old_need.best_fitness == new_need.best_fitness 
        # && (count_gates(circ, num_gates ) < count_gates(all_circuits[old_need.circuit_key], num_gates ))
      add_need = true
    end
    if add_need
      res=pop!(all_circuits,old_need.circuit_key,nothing)   # Remove old key so that former circuit won't be used in future circuits
      println("pop! result: ",res)
      new_key = "circ_"*string(new_circuit_key( key_array ))
      println("new_key: ",new_key)
      circ.key = new_key
      new_need.circuit_key = new_key
      needs[output_goal] = new_need
      all_circuits[new_key] = circ
      println("new need: ",needs[output_goal] )
      println("new circ: ",all_circuits[new_key])
    end
  end
end        

# Runs if function comp_experiment() is used to construct all_circuits
# But does not really give the desired functionality.
function random_comp_circuit( all_circuits::Dict{String,CC}, numinputs::Int64, ncircuits::Int64,
    funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs = default_funcs(numinputs)
  end
  new_circuits = CC[]
  for i = 1:ncircuits
    new_circ = deepcopy(all_circuits[rand(keys(all_circuits))])
    len = length(new_circ.inputs)
    #println("len: ",len,"  rc new_circ: ",new_circ)
    for j = 1:len
      new_circ.inputs[j] = rand(1:len+i-1)
    end
    push!(new_circuits,new_circ)
  end
  #println("rc new_circuits: ",new_circuits)
  CC(new_circuits,RandomPermutation(numinputs).data,1)
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
  new_circ = CGP.CC(gate_circuits,RandomPermutation(numinputs).data,1)
end

function check_need!( needs::Dict{Array{UInt16,1},Need}, all_circuits::Vector{CGP.CC}, circ::CGP.CC )
  #circ = all_circuits[circ_index]
  println("check need! circ: ",circ)
  @assert circ.numoutputs == 1
  numinputs = length(circ.inputs)
  perms = permutations(circ.inputs)
  for p in perms
    pcirc = deepcopy(circ)
    pcirc.inputs = p
    println("calling execute! from check_need")
    result = execute!(pcirc,construct_context(numinputs),recurse_step,needs,recurse_limit=recurse_limit)
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

function count_gates( circ::Cg, num_gates::Vector{Int64} )
  return 1
  num_gates[1] += 1
  #println("Cg cg: ",num_gates[1])
  return num_gates[1]
end

function count_gates( circ::Cc, num_gates::Vector{Int64} )
  return 1
  num_gates[1] += 1
  println("CC cg: ",num_gates[1],"  length(circ.circuits): ",length(circ.circuits))
  i = 1
  while i <= length(circ.circuits) && num_gates[1] < 20
    println("i: ",i)
    count_gates( circ.circuits[i], num_gates )
    i += 1
  end
  return num_gates[1]
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

function active_list( cc::Chromosome )
  p = cc.params
  output_values(cc)  # Execute chromosome to set nodes to active
  active_lst = zeros(Int64,p.numinputs+p.numinteriors)
  for i = 1:(p.numinputs+p.numinteriors)
    active_lst[i] = cc[i].active ? 1 : 0
  end
  active_lst
end

function active_list_random( p::Parameters, funcs::Vector{Func}, nreps::Int64 )
  active_lst = zeros(Int64,p.numinputs+p.numinteriors)
  for i = 1:nreps
    cc = random_chromosome( p, funcs )
    active_lst .+= active_list(cc)
  end
  active_lst
end

function active_list_random_df( p::Parameters, funcs::Vector{Func}, nreps::Int64 )
  active_list = active_list_random( p, funcs, nreps )
  inactive_list = nreps .- active_list
  df = DataFrame( :node=>collect(1:(p.numinputs+p.numinteriors)), :active=>active_list, :inactive=>inactive_list )
  df
end

# Convert a chromsome to a CompCircuit
# The inputs of the CompCircuit are always [1,2....]
function chromosome_to_circuit( ch::Chromosome )
  inputs = Int64[ inp.index for inp in ch.inputs ]
  gates = CGP.CC[ CGP.Cg( i.func, Vector{Int64}(i.inputs), 1 ) for i in ch.interiors ]
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
