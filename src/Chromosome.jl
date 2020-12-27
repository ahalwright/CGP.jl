using DataFrames
using CSV
import Base.getindex
export Chromosome, print_chromosome, getindex, random_chromosome, mutate_chromosome!, mutate_all
export num_mutate_locations, set_active_to_false, fraction_active, check_recursive, node_values
export output_values, number_active, number_active_gates, hamming_distance, hamming, deactivate_chromosome!
export copy_chromosome!, mutational_robustness, fault_tolerance_fitness
export build_chromosome, Input_node, Int_node, Output_node, print_build_chromosome, circuit_code, circuit_int
export circuit, print_circuit
export circuit_distance, remove_inactive, count_circuits
export code_to_circuit

mutable struct Chromosome
    params::Parameters
    inputs::Vector{InputNode}
    interiors::Vector{InteriorNode}
    outputs::Vector{OutputNode}
    fitness::Real
    robustness::Real
end

CPopulation = Vector{Chromosome}

function Chromosome(p::Parameters)
    inputs = Array{InputNode}( undef, p.numinputs)
    interiors = Array{InteriorNode}( undef, p.numinteriors )
    outputs = Array{OutputNode}( undef, p.numoutputs)
    fitness = 0.0
    robustness = 0.0
    return Chromosome(p, inputs, interiors, outputs, fitness, robustness )
end

function random_chromosome(p::Parameters, funcs::Vector{Func})
    #println("Creating random Chromosome")
    c = Chromosome(p)

    for index = 1:p.numinputs
        c.inputs[index] = InputNode(index)
        c.inputs[index].active = false
        c.inputs[index].cache = 0
    end
    for index = (p.numinputs +1):(p.numinputs+p.numinteriors)
        maxindex = index - 1 
        minindex = max(1,index-p.numlevelsback)
        #println("index: ",index,"  minindex: ",minindex,"  maxindex: ",maxindex)
        func = funcs[rand(1:end)]
        inputs = Array{Int64}(undef, func.arity)
        for i = 1:func.arity
            inputs[i] = rand(minindex:maxindex)
            #println("i: ",i,"  inputs[i]: ",inputs[i])
        end
        c.interiors[index - p.numinputs] = InteriorNode(func, inputs)
    end
    minindex = max(1,p.numinputs + p.numinteriors - p.numlevelsback + 1)
    maxindex = p.numinputs + p.numinteriors
    #println("output (minindex, maxindex): ",(minindex,maxindex))
    for i = 1:length(c.outputs)
        #index = rand(minindex:maxindex)
        index = p.numinputs + p.numinteriors + i - p.numoutputs   # use the last numoutputs interiors
        #println("output ",i," index: ",index)
        c.outputs[i] = OutputNode(index)
        c[index].active = true  #  Output nodes are always active
    end
    set_active_to_false(c)
    return c
end

function set_active_to_false( c::Chromosome )
  for index = 1:(c.params.numinputs+c.params.numinteriors)
    c[index].active = false   
    c[index].cache = MyInt(0)
  end
end

# mutates chromosome c by changing a random node or function
# if length(funcs) == 1, then only connections are mutated
# Returns new chromosome c and Bool active which is true if the node mutated was active
function mutate_chromosome!( c::Chromosome, funcs::Vector{Func}, mutate_location::Int64=0 )
  num_mutate_locs = num_mutate_locations( c, funcs )
  if mutate_location == 0 
    mutate_location = rand(1:num_mutate_locs)
  end
  #println("num_mutate_locations: ",num_mutate_locs,"  mutate_location: ",mutate_location)
  interiors_inputs_list = [ c.interiors[i].inputs for i = 1:c.params.numinteriors ]
  num_inputs_list = map(length, interiors_inputs_list)  # A list of the number of inputs of the interior node functions
  #println("num_inputs_list: ",num_inputs_list)  
  num_funcs_to_mutate = length(funcs) > 1 ? c.params.numinteriors : 0
  if mutate_location <= num_funcs_to_mutate   # mutate a func
    #println("mutate a func ml:",mutate_location)
    active = c[mutate_location].active
    new_func_index = rand(1:length(funcs))
    #println("new_funct_index: ",new_func_index,"  funcs[new_func_index].func: ",funcs[new_func_index].func)
    #println("  c.interiors[mutate_location].func: ",c.interiors[mutate_location].func.func)
    #println("condition: ", funcs[new_func_index].func == c.interiors[mutate_location].func.func )
    while( funcs[new_func_index].func == c.interiors[mutate_location].func.func )
      new_func_index = rand(1:length(funcs))
      #println("while new_funct_index: ",new_func_index)
    end
    #println("new_func_index: ",new_func_index)
    old_func_arity = c.interiors[mutate_location].func.arity
    old_inputs = c.interiors[mutate_location].inputs
    new_func = funcs[new_func_index]
    new_inputs = Array{Int64}(undef, new_func.arity)
    for i = 1:new_func.arity
      if i <= length(old_inputs)
        new_inputs[i] = old_inputs[i]
      else
        maxindex = mutate_location + c.params.numinputs - 1
        minindex = max(1,mutate_location + c.params.numinputs - c.params.numlevelsback)        
        #println("mut_location: ",mutate_location,"  minindex: ",minindex,"  maxindex: ",maxindex) 
        new_inputs[i] = rand(minindex:maxindex)
      end
    end 
    c.interiors[mutate_location] = InteriorNode(new_func,new_inputs)
    #println("new interior node: ",c.interiors[mutate_location])
    #println("new interior node: ",c[mutate_location+c.params.numinputs])
  elseif mutate_location <= num_funcs_to_mutate + sum(num_inputs_list)  # mutate an interior
    #println("mutate a interior node ml:",mutate_location)
    interior_mutate_location = mutate_location - num_funcs_to_mutate
    #println("interior mutate_location: ",interior_mutate_location)
    i = 1  # i will be the index of the interior node one of whose inputs will be mutated
    num_inputs_sum = num_inputs_list[1]
    #println("i: ",i,"  num_inputs_sum: ",num_inputs_sum)
    while(num_inputs_sum < interior_mutate_location)
      i += 1
      num_inputs_sum += num_inputs_list[i]
      #println("i: ",i,"  num_inputs_sum: ",num_inputs_sum)
    end
    #println("i: ",i,"  c.interiors[i]: ",c.interiors[i])
    active = c.interiors[i].active
    j  = num_inputs_sum - interior_mutate_location + 1  # j is the index of the input to mutate
    #println("(i,j): ",(i,j),"  num_inputs_sum: ",num_inputs_sum)
    maxindex = i + c.params.numinputs -1
    minindex = max(1, i + c.params.numinputs - c.params.numlevelsback)
    intindex = rand(minindex:maxindex)
    if maxindex > minindex
      while( intindex == c.interiors[i].inputs[j] )
        intindex = rand(minindex:maxindex)
      end
    end 
    #println("minindex: ",minindex,"  maxindex: ",maxindex,"  intindex: ",intindex)
    c.interiors[i].inputs[j] = intindex
    #println("(i,j): ",(i,j),"  c.interiors[i].inputs: ",c.interiors[i].inputs)
  elseif mutate_location > num_funcs_to_mutate + sum(num_inputs_list)   # mutate an output
    # Mutating an output node should never happen because random chromosomes always use that last interior nodes as input,
    #   and num_mutate_locations() returns a value that does not include mutation of output nodes.
    error("mutate_chromosome!() is mutating an output node which should never happen")
    #println("mutate the output: ",mutate_location - num_funcs_to_mutate - sum(num_inputs_list ))
    active = true   # outputs are always active
    minindex = max(1, num_funcs_to_mutate + c.params.numinputs - c.params.numlevelsback + 1 ) 
    maxindex = num_funcs_to_mutate + c.params.numinputs
    if maxindex == mutate_location
      print_chromosome(c)
      exit()
    end
    #println("minindex: ",minindex,"  maxindex: ",maxindex)
    if maxindex > minindex   # no change if there is only one alternative for the output node index
      outindex = rand(minindex:maxindex)
      #println("outindex: ",outindex)
      #println("num_funcs_to_mutate: ", num_funcs_to_mutate )
      #println("sum(num_inputs_list): ", sum(num_inputs_list))
      while outindex == c.outputs[mutate_location - num_funcs_to_mutate - sum(num_inputs_list)].input
        outindex = rand(minindex:maxindex)
      end
      c.outputs[mutate_location - num_funcs_to_mutate - sum(num_inputs_list)].input = outindex
    end
  end
  #println("active: ",active)
  set_active_to_false(c)
  (c,active)
end

# Mutates c in all possible ways and returns a list of the outputs of the mutated chromosomes
#    or the pair (avg_robustness, evolvability)
# Deterministic if all functions in default_funcs() have the same arity
# If robustness_only==true returns the pair: (avg_robustness, evolvability)
# Otherwise returns vector of either outputs, or chromosomes, or (outputs, chromosomes) pairs
function mutate_all( c::Chromosome, funcs::Vector{Func}; 
      robustness_only::Bool=false, output_outputs::Bool=true, output_chromosomes::Bool=false )
  #println("mutate_all: numlevelsback: ", c.params.numlevelsback )
  sav_c = deepcopy(c)
  if robustness_only
    output_outputs = output_chromosomes = false
    robustness_sum = 0.0; robustness_count = 0
    result = Vector{MyInt}[]   # Save outputs to give a measure of evolvability
  elseif output_outputs && !output_chromosomes
    result = Vector{MyInt}[]
  elseif !output_outputs && output_chromosomes 
    result = Chromosome[]
  elseif output_outputs && output_chromosomes 
    result = Tuple{Vector{MyInt},Chromosome}[]
  end
  context = construct_context(c.params.numinputs)
  orig_output = output_values(c)
  num_mutate_locs = num_mutate_locations( c, funcs )
  interiors_inputs_list = [ c.interiors[i].inputs for i = 1:c.params.numinteriors ]
  num_inputs_list = map(length, interiors_inputs_list)
  num_funcs_to_mutate = length(funcs) > 1 ? c.params.numinteriors : 0
  for mutate_location = 1:num_mutate_locs
    new_c = deepcopy(c)
    if mutate_location <= num_funcs_to_mutate   # mutate a func
      #println("mutate a func ml:",mutate_location,"  node: ",c.interiors[mutate_location])
      active = c[mutate_location].active
      for new_func_index = 1:length(funcs)
        # Does not generate all alternative mutations in the very unusual case when the new function has greater arity
        #    than the old function
        #println("new_func: ",funcs[new_func_index].func,"  old_func: ", c.interiors[mutate_location].func.func )
        if funcs[new_func_index].func == c.interiors[mutate_location].func.func
          continue   # skip
        end
        #println("new_funct_index: ",new_func_index,"  funcs[new_func_index].func: ",funcs[new_func_index].func)
        #println("  c.interiors[mutate_location].func: ",c.interiors[mutate_location].func.func)
        #println("condition: ", funcs[new_func_index].func == c.interiors[mutate_location].func.func )
        #println("new_func_index: ",new_func_index)
        old_func_arity = c.interiors[mutate_location].func.arity
        old_inputs = c.interiors[mutate_location].inputs
        new_func = funcs[new_func_index]
        new_inputs = Array{Int64}(undef, new_func.arity)
        for i = 1:new_func.arity
          if i <= length(old_inputs)
            new_inputs[i] = old_inputs[i]
          else
            maxindex = mutate_location + c.params.numinputs - 1
            minindex = max(1,mutate_location + c.params.numinputs - c.params.numlevelsback)        
            #println("mut_location: ",mutate_location,"  minindex: ",minindex,"  maxindex: ",maxindex) 
            new_inputs[i] = rand(minindex:maxindex)
          end
        end 
        new_c = deepcopy(c)
        new_c.interiors[mutate_location] = InteriorNode(new_func,new_inputs)
        deactivate_chromosome!(new_c)
        new_output = execute_chromosome(new_c,context)
        if robustness_only
          robustness_sum = (new_output == orig_output) ? robustness_sum + 1 : robustness_sum
          robustness_count+=1
          push!(result,new_output)
        elseif output_outputs && !output_chromosomes
          push!(result,new_output)
        elseif !output_outputs && output_chromosomes
          push!(result,new_c)
        elseif output_outputs && output_chromosomes
          push!(result,(new_output,new_c))
        end
      end
      #println("new interior node: ",new_c.interiors[mutate_location])
    elseif mutate_location <= num_funcs_to_mutate + sum(num_inputs_list)  # mutate an interior
      #println("mutate an interior node ml:",mutate_location, "  node: ",c.interiors[mutate_location-num_funcs_to_mutate])
      interior_mutate_location = mutate_location - num_funcs_to_mutate
      #println("interior mutate_location: ",interior_mutate_location) 
      i = 1  # i will be the index of the interior node one of whose inputs will be mutated
      num_inputs_sum = num_inputs_list[1]
      #println("i: ",i,"  num_inputs_sum: ",num_inputs_sum)
      while(num_inputs_sum < interior_mutate_location)
        i += 1
        num_inputs_sum += num_inputs_list[i]
        #println("i: ",i,"  num_inputs_sum: ",num_inputs_sum)
      end
      #println("i: ",i,"  c.interiors[i]: ", c.interiors[i])
      j  = num_inputs_sum - interior_mutate_location + 1  # j is the index of the input to mutate
      maxindex = i + c.params.numinputs -1
      minindex = max(1, i + c.params.numinputs - c.params.numlevelsback)
      #println("(i,j): ",(i,j),"  num_inputs_sum: ",num_inputs_sum,"  minidex: ",minindex,"  maxindex: ",maxindex)
      for intindex = minindex:maxindex
        #println("intindex: ",intindex,"  c.interiors[i].inputs[j]: ",c.interiors[i].inputs[j])
        if intindex == c.interiors[i].inputs[j]
          #println("continue")
          continue
        end
        new_c = deepcopy(c)
        new_c.interiors[i].inputs[j] = intindex
        #println("intindex: ",intindex,"  (i,j): ",(i,j),"  new_c.interiors[i].inputs: ",new_c.interiors[i].inputs)
        deactivate_chromosome!(new_c)
        new_output = execute_chromosome(new_c,context)
        if robustness_only
          robustness_sum = (new_output == orig_output) ? robustness_sum + 1 : robustness_sum
          robustness_count+=1
          push!(result,new_output)
        elseif output_outputs && !output_chromosomes
          push!(result,new_output)
        elseif !output_outputs && output_chromosomes
          push!(result,new_c)
        elseif output_outputs && output_chromosomes
          push!(result,(new_output,new_c))
        end
      end
    elseif mutate_location > num_funcs_to_mutate + sum(num_inputs_list)   # mutate an output
      # Mutating an output node should never happen because random chromosomes always use that last interior nodes as input,
      #   and num_mutate_locations() returns a value that does not include mutation of output nodes.
      error("mutate_all() is mutating an output node which should never happen")
      #println("mutate the output: ",mutate_location - num_funcs_to_mutate - sum(num_inputs_list ))
      #active = true   # outputs are always active
      minindex = max(1, num_funcs_to_mutate + c.params.numinputs - c.params.numlevelsback + 1 ) 
      maxindex = num_funcs_to_mutate + c.params.numinputs
      if maxindex == mutate_location
        print_chromosome(c)
        exit()
      end
      #println("minindex: ",minindex,"  maxindex: ",maxindex)
      if maxindex > minindex   # no change if there is only one alternative for the output node index
        for outindex = minindex:maxindex
          #outindex = rand(minindex:maxindex)
          if outindex == c.outputs[mutate_location - num_funcs_to_mutate - sum(num_inputs_list)].input 
            continue
          end
        end
        println("outindex: ",outindex)
        new_c.outputs[mutate_location - num_funcs_to_mutate - sum(num_inputs_list)].input = outindex
        new_output = execute_chromosome(new_c,context)
        push!(result,new_c)
      end
    end
  end
  if robustness_only
    (robustness_sum/robustness_count, length(unique(result))/length(result)) # pair of avg robustenss and evolvability
  else
    result
  end
end

function num_mutate_locations( c::Chromosome, funcs::Vector{Func} )
  num_funcs_to_mutate = length(funcs) > 1 ? c.params.numinteriors : 0
  interiors_inputs_list = [ c.interiors[i].inputs for i = 1:c.params.numinteriors ]
  outputs_input_list = [ c.outputs[i].input for i = 1:c.params.numoutputs ]
  #println("int_list: ", interiors_inputs_list, "   out_list: ",outputs_input_list)
  num_inputs_list = map(length, interiors_inputs_list)
  #num_inputs = sum(num_inputs_list) + c.params.numoutputs
  num_inputs = sum(num_inputs_list)
  #println("num_inputs_list: ",num_inputs_list,"  num_inputs: ",num_inputs)
  if length(funcs) == 1
    num_mutate_locations = num_inputs 
  else
    num_mutate_locations = num_funcs_to_mutate + num_inputs 
  end
  num_mutate_locations
end

function print_chromosome(f::IO, c::Chromosome )
  for i = 1:(c.params.numinputs + c.params.numinteriors + c.params.numoutputs)
    println(f,"i: ",i,"  ",c[i])
  end
  check_recursive(c)
end

function print_chromosome( c::Chromosome )
  for i = 1:(c.params.numinputs + c.params.numinteriors + c.params.numoutputs)
    println("i: ",i,"  ",c[i])
  end
  check_recursive(c)
end

function check_recursive( c::Chromosome )
  for i = (c.params.numinputs + 1) : c.params.numinteriors
    for j = 1:length(c[i].inputs)
      #println("i: ",i,"  inputs[",j,"] = ",c[i].inputs[j])
      if i == j 
        println("recursive ")
        exit()
      end
    end
  end
end

function number_active_old( c::Chromosome )
  num_active = 0
  for i = 1:c.params.numinputs
    if c.inputs[i].active
      num_active += 1
    end
  end
  for i = 1:c.params.numinteriors
    if c.interiors[i].active
      num_active += 1
    end
  end
  num_active
end

function number_active( c::Chromosome )
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
    context = construct_context(c.params.numinputs)
    execute_chromosome(c,context)
  end
  num_act_in = reduce(+,[c.inputs[i].active for i = 1:length(c.inputs)])
  num_act_int =  reduce(+,[c.interiors[i].active for i = 1:length(c.interiors)])
  num_act_in + num_act_int
end

function number_active_gates( c::Chromosome )
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
    context = construct_context(c.params.numinputs)
    execute_chromosome(c,context)
  end
  num_act_int =  reduce(+,[c.interiors[i].active for i = 1:length(c.interiors)])
end

function node_values( c::Chromosome )
  output_values = [c[c.outputs[i].input].cache for i = 1:length(c.outputs)]
  input_values = [c.inputs[i].cache for i = 1:length(c.inputs)]
  interior_values = [c.interiors[i].cache for i = 1:(length(c.interiors))]
  #output_values = [c.interiors[i].cache for i = (length(c.interiors)-length(c.outputs)+1):length(c.interiors)]
  (input_values,interior_values,output_values)
end

function output_values( c::Chromosome )
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
    context = construct_context(c.params.numinputs)
    return execute_chromosome(c,context)
  else
    return [c[c.outputs[i].input].cache for i = 1:length(c.outputs)]
  end
end

function deactivate_chromosome!( c::Chromosome )
  for in in c.inputs
    in.active = false
  end
  for int in c.interiors
    int.active = false
  end
  c
end

function fraction_active( c::Chromosome )
  number_active(c)/(c.params.numinputs+c.params.numinteriors)
end

# Returns a chromosome with the same output values as c and with all 
#  inactive gate nodes removed.
function remove_inactive( c::Chromosome )
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
    context = construct_context(c.params.numinputs)
    execute_chromosome(c,context)
  end
  new_ints = InteriorNode[]
  numinputs = length(c.inputs)
  numints = length(c.interiors)
  int_inds = collect(1:numints)
  # if i refers to a node in c, then map_indices[i] references 
  #   the corresponding node in the new chromosome
  map_indices = zeros(Int64, numinputs+numints )
  active_gate_indices = int_inds[ [c.interiors[i].active for i = 1:numints]]
  #map_indices = ttt(numinputs,numints,active_gate_indices)
  map_indices = collect(1:numinputs)
  j = 1
  for i = 1:numints
    if active_gate_indices[j] > i
      push!(map_indices,0)
    else
      push!(map_indices,j+numinputs)
      j+= 1
    end
  end
  #println("map_indices: ",map_indices)
  for j in active_gate_indices
    new_inputs = map(x->map_indices[x],c.interiors[j].inputs)
    #println("j: ",j,"  new_inputs: ",new_inputs)
    new_int = InteriorNode( c.interiors[j].func, new_inputs )
    push!(new_ints,new_int)
  end
  new_outs = [ OutputNode( map_indices[c.outputs[i].input])  for i = 1:length(c.outputs) ]
  p = c.params
  new_p = Parameters( p.numinputs, p.numoutputs, length(new_ints), p.numlevelsback )
  nc = Chromosome( new_p, c.inputs, new_ints, new_outs, c.fitness, c.robustness )
  @assert output_values(nc) == output_values(c)
  nc
end             

function ttt(numinputs,numints,active_gate_indices)
  map_indices = collect(1:numinputs)
  j = 1
  for i = 1:numints
    if active_gate_indices[j] > i
      push!(map_indices,0)
    else
      push!(map_indices,j+numinputs)
      j+= 1
    end
  end
  map_indices
end

function getindex(c::Chromosome, index::Integer)
    if index <= c.params.numinputs   # input node
        return c.inputs[index]
    end

    if index > c.params.numinteriors + c.params.numinputs  # output node
        return c.outputs[index - c.params.numinteriors - c.params.numinputs]
    end

    return c.interiors[index-c.params.numinputs]  # interior node
end

# Number of bits where x and y differ.  x and y are interpreted as bit strings
function hamming( x::MyInt, y::MyInt )
  xr = xor( x, y )
  result = 0
  my_one = convert(MyInt,1)
  my_zero = convert(MyInt,0)
  while xr != my_zero
    result +=  xr & my_one
    xr >>= 1
  end
  result
end

# Hamming distance between MyInts which is between 0 and 1.
# If hamming(x,y) >= 2^numinputs/2, then hamming_distance(x,y,numinputs) = 1.0
function hamming_distance( x::MyInt, y::MyInt, numinputs::Int64 )
  result = hamming(x,y)/2^numinputs
end

# Hamming distance between Goals which is between 0 and 1.
function hamming_distance( x::Vector{MyInt}, y::Vector{MyInt}, numinputs::Int64 )
  @assert length( x ) == length( y )
  sum( hamming_distance( x[i], y[i], numinputs ) for i = 1:length(x))/length(x)
end

# Change the fields of chromosome c to be the fields of chromosom c_to_copy
function copy_chromosome!( c::Chromosome, c_to_copy::Chromosome )
  c.params = c_to_copy.params
  c.inputs = c_to_copy.inputs
  c.interiors = c_to_copy.interiors
  c.outputs = c_to_copy.outputs
  c.interiors = c_to_copy.interiors
end

# Computes the fraction of mutations that do not change the output
# Not really accurate because it counts one mutation per num_mutate_location(), while in fact
#    there may be more than one random choice per num_mutate_location()
# Superceded by calling   mutate_all( c, funcs, robustness_only=true ) which returns the pair (robustness,evolvability)
# If active_only==true, then only mutations that change active nodes are considered
# If active_only==false, then all mutations are considered
# Note that results are partially random since the result of mutations at a specific location can be different.
function mutational_robustness( c::Chromosome, funcs::Vector{Func}; active_only::Bool=false ) 
  context = construct_context( c.params.numinputs )
  prev_out = execute_chromosome( c, context )  
  count_no_change = 0
  count_no_change = 0.0
  num_mut_locs = 0
  for j = 1:num_mutate_locations( c, funcs )
    cc = deepcopy(c)
    (cc,active) = mutate_chromosome!( cc, funcs, j )
    out_c = execute_chromosome( cc, context )
    #println("j: ",j,"  active: ",active)
    #print_chromosome( cc )
    if !active_only || active # consider all mutations if !active_only, only active mutations if active_only
      #out_c = execute_chromosome( cc, context )
      count_no_change += ((out_c == prev_out) ? 1 : 0)
      if active_only && active
        num_mut_locs += 1
      end
    end
    #println("j: ",j,"  active: ",active,"  no_change: ",(out_c == prev_out),"  count_no_change: ",count_no_change)
  end
  #println("num_mut_locs: ",num_mut_locs,"  num_mutate_locations(): ",num_mutate_locations(c,funcs),"  count_no_change: ",count_no_change)
  #@assert num_mut_locs == num_mutate_locations(c,funcs)
  if !active_only
    num_mut_locs = num_mutate_locations( c, funcs )
  end
  count_no_change/num_mut_locs
end

# Implements equation 3.3 of Macia and Sole (2009)
#  The average Hamming deviation of the output from the unperturbed output under perturbation of the output of each node
function fault_tolerance_fitness( c::Chromosome )
  numinputs = c.params.numinputs
  numoutputs = c.params.numoutputs
  numints = c.params.numinteriors
  #println("(numinputs,numoutputs,numints):",(numinputs,numoutputs,numints))
  context = construct_context( c.params.numinputs )
  outputs = output_values( c )
  #println("outputs: ",outputs)
  distances = fill(0.0,numints)
  for i = (1+numinputs):(numinputs+numints)    
    out_ft = execute_chromosome_ft(c,context,i)
    #println("i: ",i,"  out_ft: ",out_ft)
    #println("hamming dis: ",hamming_distance(execute_chromosome_ft(c,context,i)[1], outputs[1], numinputs))
    distances[i-numinputs] =  
      sum( hamming_distance(execute_chromosome_ft(c,context,i)[k], outputs[k], numinputs) for k = 1:numoutputs)
  end
  #println("distances: ",distances)
  1.0-sum(distances)/numints/numoutputs
end

#Assumes nodearity==2
mutable struct Gate_node
  node_index::Int64
  func::Func
  input1::Int64
  input2::Int64
end

# Creates a Chromosome (circuit) from the concise format.
# Example:  
# julia>circuit((1,2,3), ((4,OR,1,2), (5,AND,2,3), (6,XOR,4,5)))
#    creates a ciruit with 3 input nodes and 3 gate nodes and 1 implicit output node.
# gate nodes have 4 fields: 
#    node_index:  should be the index of the node.  Not used in construcing the circuit, but checked for correctness
#    node_function:  See Func.jl for options
#    node_input1:  an integer in the range from 1 to the index of the node minus 1
#    node_input2:  an integer in the range from 1 to the index of the node minus 1
# Assumes nodearity==2 and numoutputs == 1
function circuit( inputs::Tuple, gates::Tuple;
    levsback::Int64=length(inputs)+length(gates))
  numinputs = length(inputs)
  errors = false
  for i = 1:length(gates)
    if gates[i][1] != i+numinputs
      println("index of gate ",i," not correct: it should be: ",i+numinputs)
      errors = true
    end
  end
  for i = 1:length(gates)
    if !((1 <= gates[i][3])  && (gates[i][3] < numinputs+i))
      error("illegal first gate input for gate: ",i)
    end
  end
  for i = 1:length(gates)
    if !((1 <= gates[i][4])  && (gates[i][4] < numinputs+i))
      error("illegal second gate input for gate: ",i)
    end
  end
  p = Parameters( numinputs=length(inputs), numoutputs=1, numinteriors=length(gates), numlevelsback=levsback )
  in_nodes = [InputNode(i) for i in inputs]   
  gate_nodes = [InteriorNode(g[2], [g[3],g[4]] ) for g in gates ]
  out_nodes = [OutputNode( length(inputs)+length(gates ) ) ]
  c = Chromosome( p, in_nodes, gate_nodes, out_nodes, 0.0, 0.0 )
  if errors
    println("correct input to this function should have been: ")
    print_circuit(c)
  end
  c
end
  
# Outputs the concise format of Chromosome (circuit) c to IO stream f
# print_node_tuple() and gate_tuple() are defined in Node.jl
function print_circuit( f::IO, c::Chromosome; include_fitness::Bool=false )  
  @assert c.params.nodearity==2
  @assert c.params.numoutputs==1
  #@assert c.params.numinputs>1
  print(f,"circuit(")
  print_node_tuple(f, c.inputs )
  print(f,",")
  gate_tuple(f, c.interiors, c.params.numinputs )
  if !include_fitness
    println(f,")")
  else
    println(f,", 0.0)")
  end
  #=
  if typeof(f) == IOStream  # Don't close Base.stdout since this kills julia
    close(f)
  end
  =#
end

function print_circuit( c::Chromosome; include_fitness::Bool=false )
  print_circuit( Base.stdout, c, include_fitness=include_fitness )
end    

# The following struct defintions and function definitions build_chromosome() and print_build_chromosome()
#   have mostly been replaced by the above definitions of circuit() and print_circuit()
mutable struct Int_node
  func::Func
  inputs::Vector{Int64}
end

mutable struct Input_node
  index::Int64
end

mutable struct Output_node
  input::Integer
end

# Example call: build_chromosome( [Input_node(1),Input_node(2)], [Int_node(OR,[1,2]),Int_node(AND,[2,3])],[Output_node(4)])
# Depreciated 12/18/20
function build_chromosome( input_nodes::Vector{Input_node}, interior_nodes::Vector{Int_node}, output_nodes::Vector{Output_node};
    levsback::Int64=length(input_nodes)+length(interior_nodes))
  num_in = length(input_nodes)
  num_ints = length(interior_nodes)
  num_outs = length(output_nodes)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=levsback )
  in_nodes = [InputNode(in_node.index) for in_node in input_nodes]
  int_nodes = [InteriorNode(int_node.func, int_node.inputs) for int_node in interior_nodes]
  out_nodes = [OutputNode(outnode.input) for outnode in output_nodes]
  Chromosome( p, in_nodes, int_nodes, out_nodes, 0.0, 0.0 )
end

# Example calls:  build_chromosome((1,2), ((OR,[1,2]),(AND,[2,3])),(4,))
# Example calls:  build_chromosome((1,2), ((OR,[1,2]),(AND,[2,3])),(4,),0.0)
# Depreciated 12/18/20
function build_chromosome( inputs::Tuple, ints::Tuple, outs::Tuple, fitness::Float64=0.0; 
    levsback::Int64=length(inputs)+length(ints))
  num_in = length(inputs)
  num_ints = length(ints)
  num_outs = length(outs)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=levsback )
  in_nodes = [InputNode(in_index) for in_index in inputs]
  int_nodes = [InteriorNode(int_pair[1], int_pair[2]) for int_pair in ints]
  out_nodes = [OutputNode(out_index) for out_index in outs]
  Chromosome( p, in_nodes, int_nodes, out_nodes, fitness, 0.0 )
end

function print_build_chromosome( f::IO, c::Chromosome; include_fitness::Bool=false )
  print(f, "build_chromosome(")
  print_node_tuple(f, c.inputs )
  print(f,",")
  print_node_tuple(f, c.interiors )
  print(f,",")
  print_node_tuple(f, c.outputs )
  if include_fitness
    println(f,", ",c.fitness)
  end
  println(f,")")
  if typeof(f) == IOStream  # Don't close Base.stdout since this kills julia
    close(f)
  end
end

function print_build_chromosome( c::Chromosome )
  print_build_chromosome( Base.stdout, c)
end          

# Returns a vector of integers which is a characterization of the chromosome for the parameters c.params
# The length of the result should be p.numinteriors*(1+p.nodearity).
# Not correct in that codes for gate inputs are sometimes too large.
function circuit_code( c::Chromosome )
  result = Int64[]
  funcs = default_funcs(c.params.numinputs)
  if length(c.interiors) == 0
    return Int64[]
  end
  for i = 1:c.params.numinteriors
    # Determine index of gate function
    j = 0
    while j < length(funcs) && funcs[j+1].func != c.interiors[i].func.func
      j += 1
    end
    if j > length(funcs)
      error("illegal function in function circuit_code()")
    end
    #println("i: ",i,"  function code: ",j)
    push!(result,j)
    # Add interior node inputs to result
    index = i + c.params.numinputs  # node index in chromosome
    for k = 1:length(c.interiors[i].inputs)
      code = max(c.interiors[i].inputs[k]-(index-c.params.numlevelsback),0)
      println("i: ",i,"  k: ",k,"  code: ",code)
      push!(result, code)
    end
  end
  # Note:  assumes that the input fields of output nodes are set by default and thus don't need to be recorded
  result
end

# Not correct
# Converts a circuit code to a Chromosome.   Inverse function to circuit_code().
function code_to_circuit( code::Vector{Int64}, p::Parameters )
  funcs = default_funcs( p.numinputs )
  interior_nodes = InteriorNode[]
  for i = 1:p.numinteriors
    j = 3*(i-1)+1   # index of code component
    func = funcs[code[j] + 1]
    index = i + p.numinputs  # node index in Chromosome
    input1 = code[j+1] + max(index-p.numlevelsback,0)
    input2 = code[j+2] + max(index-p.numlevelsback,0)
    #println("i: ",i,"  j: ",j,"  func: ",func,"  input1: ",input1,"  input2: ",input2)
    push!( interior_nodes, InteriorNode(func,[input1,input2]) )
  end
  input_nodes = [ InputNode(j) for j = 1:p.numinputs ]
  output_nodes = [ OutputNode(j) for j = (p.numinputs+p.numinteriors-p.numoutputs+1):(p.numinputs+p.numinteriors) ]
  c = Chromosome(p, input_nodes, interior_nodes, output_nodes, 0.0, 0.0 )
end

# An integer that characterizes the chromosome
# Doesn't overflow for 11 gates, 8 numlevelsback
function circuit_int( c::Chromosome )
  circuit_int( circuit_code(c), c.params )
end

# Not correct  See diary12_25.txt for example.
function circuit_int( c_code::Vector{Int64}, p::Parameters )
  result = Int128(0)
  multiplier = 1
  funcs = default_funcs( p.numinputs )
  k = 1
  for i = 1:p.numinteriors
    multiplier = length(funcs)
    result = c_code[k] + result*multiplier
    print("i: ",i,"  k: ",k,"  multiplier: ",multiplier,"  c_code[k]: ",c_code[k])
    println("  result: ",result)
    k += 1
    multiplier = min(p.numlevelsback,i-1+p.numinputs)
    for j = 1:p.nodearity
      #println("result: ",result,"  multiplier: ",multiplier,"  c_code[k]: ",c_code[k])
      result = c_code[k] + result*multiplier
      print("i: ",i,"  j: ",j,"  k: ",k,"  multiplier: ",multiplier,"  c_code[k]: ",c_code[k])
      println("  result: ",result)
      k += 1
    end
  end
  result
end
function int_to_circuit_code( c_int::Integer, p::Parameters )
  c_int = Int128(c_int)
  c_code = zeros(Int64,3*p.numinteriors)
  k = 3*p.numinteriors
  for i = p.numinteriors:-1:1
    multiplier = min(p.numlevelsback,i-1+p.numinputs)
    println("i: ",i,"  multiplier: ",multiplier)
    for j = p.nodearity:-1:1
      c_int_mod= c_int % multiplier 
      #prev_c_int =c_int
      c_int ÷= multiplier
      #c_int_mod = prev_c_int - c_int
      println("i: ",i,"  j: ",j,"  c_int: ",c_int,"  c_int_mod: ",c_int_mod)
      c_code[k] = c_int_mod
      k -= 1
    end
    multiplier = length(funcs)
    println("i: ",i,"  multiplier: ",multiplier)
    c_int_mod= c_int % multiplier 
    c_int ÷= multiplier
    println("i: ",i,"  c_int: ",c_int,"  c_int_mod: ",c_int_mod)
    c_code[k] = c_int_mod
    k -= 1
  end
  c_code
end

# Hamming distance between chromosomes normalized to be between 0.0 and 1.0
function circuit_distance( c1::Chromosome, c2::Chromosome )
  code1 = circuit_code(c1)
  code2 = circuit_code(c2)
  @assert length(code1) == length(code2)
  diff_count = 0
  for i = 1:length(code1)
    diff_count += code1[i] == code2[i] ? 0 : 1
  end
  diff_count/length(code1)
end

# Return the number of circuits for parameters p and the corresponding default_funcs.
function count_circuits( p::Parameters )
  @assert p.numoutputs == 1   # Not tested for more than 1 output, but probably works in this case.
  funcs = default_funcs(p.numinputs)
  multiplier = Int128(1)
  mij = 0
  for i = 1:p.numinteriors
    mf = length(funcs)
    multiplier *= mf
    for j = 1:p.nodearity
      mij = min(p.numlevelsback,i-1+p.numinputs)
      multiplier *= mij
    end
    print("i: ",i,"  mf: ",mf,"  mij: ",mij,"  log multiplier: ",log10(multiplier))
    exp = trunc(log10(multiplier))
    fract = 10^(log10(multiplier)-exp)
    println("  exp: ",exp,"  fract: ",fract)
    @printf("  multiplier: %4.2f",fract)
    @printf("e+%2i\n",exp)
    try
      @printf("  multiplier:  %8.2e\n",multiplier)
    catch
    end
  end
  multiplier
end
