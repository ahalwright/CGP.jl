import Base.getindex
export Chromosome, print_chromosome, getindex, random_chromosome, mutate_chromosome!, mutate_all
export num_mutate_locations, set_active_to_false, fraction_active, check_recursive, node_values
export output_values, number_active, number_active_old, hamming_distance, hamming, deactivate_chromosome!
export copy_chromosome!, mutational_robustness, fault_tolerance_fitness
export build_chromosome, Input_node, Int_node, Output_node, print_build_chromosome

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
# Deterministic if all functions in default_funcs() have the same arity
# If robustness_only return a pair of avg robustness and evolvability
# Otherwise returns vector of either outputs, or chromosomes, or (outputs, chromosomes) pairs
function mutate_all( c::Chromosome, funcs::Vector{Func}; 
      robustness_only::Bool=false, output_outputs::Bool=true, output_chromosomes::Bool=false )
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

# Hamming distance which is between 0 and 1.
# If hamming(x,y) >= 2^numinputs/2, then hamming_distance(x,y,numinputs) = 1.0
function hamming_distance( x::MyInt, y::MyInt, numinputs::Int64 )
  result = hamming(x,y)/2^numinputs
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
function build_chromosome( input_nodes::Vector{Input_node}, interior_nodes::Vector{Int_node}, output_nodes::Vector{Output_node} )
  num_in = length(input_nodes)
  num_ints = length(interior_nodes)
  num_outs = length(output_nodes)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=num_ints+num_outs )
  in_nodes = [InputNode(in_node.index) for in_node in input_nodes]
  int_nodes = [InteriorNode(int_node.func, int_node.inputs) for int_node in interior_nodes]
  out_nodes = [OutputNode(outnode.input) for outnode in output_nodes]
  Chromosome( p, in_nodes, int_nodes, out_nodes, 0.0, 0.0 )
end

# Example calls:  build_chromosome((1,2), ((OR,[1,2]),(AND,[2,3])),(4,))
# Example calls:  build_chromosome((1,2), ((OR,[1,2]),(AND,[2,3])),(4,),0.0)
function build_chromosome( inputs::Tuple, ints::Tuple, outs::Tuple, fitness::Float64=0.0 )
  num_in = length(inputs)
  num_ints = length(ints)
  num_outs = length(outs)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=num_ints+num_outs )
  in_nodes = [InputNode(in_index) for in_index in inputs]
  int_nodes = [InteriorNode(int_pair[1], int_pair[2]) for int_pair in ints]
  out_nodes = [OutputNode(out_index) for out_index in outs]
  Chromosome( p, in_nodes, int_nodes, out_nodes, fitness, 0.0 )
end

function print_build_chromosome( f::IO, c::Chromosome )
  print(f, "build_chromosome(")
  print_node_tuple(f, c.inputs )
  print(f,",")
  print_node_tuple(f, c.interiors )
  print(f,",")
  print_node_tuple(f, c.outputs )
  println(f,", ",c.fitness,")")
  #print(f, ") ")
  if typeof(f) == IOStream  # Don't close Base.stdout since this kills julia
    close(f)
  end
end

function print_build_chromosome( c::Chromosome )
  print_build_chromosome( Base.stdout, c)
end          

