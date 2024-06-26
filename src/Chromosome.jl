using DataFrames
using CSV
using Printf
import Base.getindex
export PredType   # defined below to be Int64
export Chromosome, print_chromosome, getindex, random_chromosome, mutate_chromosome!, mutate_all, mutate_all_neutral
export num_mutate_locations, set_active_to_false, fraction_active, check_recursive, node_values
export output_values, number_active, number_active_gates, remove_inactive, deactivate_chromosome!
export hamming_distance, ihamming_distance, hamming 
export copy_chromosome!, mutational_robustness, fault_tolerance_fitness, number_active_old, append_chromosome
export build_chromosome, Input_node, Int_node, Output_node, print_build_chromosome, circuit_code 
export circuit, print_circuit, geno_distance
export circuit_distance, remove_inactive, count_circuits_ch, count_genotypes_ch, count_genotypes_ch_mt
export insert_gate!, delete_gate!, test_combine_complexity, combine_chromosomes, interleave_chromosomes
export enumerate_circuits_ch, chromosome_to_int, gate_int, gate_int_list, int_to_gate, int_to_chromosome
export chromosome_add_input, chromsome_add_multiplexer
export normalize_chromosome, normalize_chromosome!, count_circuits_ch_normalize, count_genotypes_table
export inputs_list
#export neutral_component_ints!
  
PredType = Int64   # Also defined in aliases.jl
#= Commented out.  The included definition is in aliases.jl.
mutable struct Chromosome
  params::Parameters
  inputs::Vector{InputNode}
  interiors::Vector{InteriorNode}
  outputs::Vector{OutputNode}
  fitness::Float64
  robustness::Union{Float64,PredType}
end
=#
CPopulation = Vector{Chromosome}

# A chromosome includes InputNodes, InteriorNodes, and OutputNodes in that order.
# This function defines an indexing on these notes.
# Thus, if Chromosome ch has 4 InputNodes, then getindex(ch,5) would return the first interior node.
#=
function getindex(c::Chromosome, i::Integer)
  if i <= c.params.numinputs   # input node
    return c.inputs[i]
  elseif i <= c.params.numinteriors 
    return c.interiors[i-c.params.numinputs]  # interior node
  else  # outputnode
    return c.outputs[i - c.params.numinteriors - c.params.numinputs]
  end
end
=#
function getindex(c::Chromosome, index::Integer)
    if index <= c.params.numinputs   # input node
        return c.inputs[index]
    end
    if index > c.params.numinteriors + c.params.numinputs  # output node
        return c.outputs[index - c.params.numinteriors - c.params.numinputs]
    end
    return c.interiors[index-c.params.numinputs]  # interior node
end

function Chromosome(p::Parameters; ident::PredType=PredType(0))
  inputs = Array{InputNode}( undef, p.numinputs)
  interiors = Array{InteriorNode}( undef, p.numinteriors )
  outputs = Array{OutputNode}( undef, p.numoutputs)
  fitness = 0.0
  if ident == PredType(0)
    robustness = 0.0
  else
    robustness = ident
  end
  return Chromosome(p, inputs, interiors, outputs, fitness, robustness )
end

function random_chromosome(p::Parameters; ident::PredType=PredType(0))
  funcs = default_funcs(p.numinputs)
  random_chromosome(p,funcs, ident=ident)
end

function random_chromosome(p::Parameters, funcs::Vector{Func}; ident::PredType=PredType(0))
  #println("Creating random Chromosome")
  c = Chromosome(p)
  for index = 1:p.numinputs
      c.inputs[index] = InputNode(index, false, MyInt(0))
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
  #println("c: ",c)
  for i = 1:length(c.outputs)
      #index = rand(minindex:maxindex)
      index = p.numinputs + p.numinteriors + i - p.numoutputs   # use the last numoutputs interiors
      #println("output i: ",i," index: ",index)
      c.outputs[i] = OutputNode(index)
      c[index].active = true  #  Output nodes are always active
  end
  set_active_to_false!(c)
  if ident != PredType(0)
    c.robustness = ident
  end
  return c
end

function set_active_to_false!( c::Chromosome )
  for index = 1:(c.params.numinputs+c.params.numinteriors)
    c[index].active = false   
    c[index].cache = MyInt(0)
  end
  c
end

# mutates chromosome c by mutating one of the gates in gate_list
function mutate_chromosome_gates!( c::Chromosome, funcs::Vector{Func}, gate_list::Vector{Int64} )
  p = c.params
  len_gts = length(gate_list)
  #func_list = deepcopy(gate_list)
  inputs_list = Int64[]
  for g in gate_list
    push!( inputs_list, g )
    push!( inputs_list, 2*g )
    push!( inputs_list, 2*g+1 )
  end
  location_list = inputs_list
  loc = rand(location_list)
  #println("mutate_chromosome_gates!(): location_list: ",location_list,"  loc: ",loc)
  mutate_chromosome!( c, funcs, loc )
end

# mutates chromosome c by changing a random node or function
# if length(funcs) == 1, then only connections are mutated
# Returns new chromosome c and Bool active which is true if the node mutated was active
# BUG!! as of 3/3/21:  Does not modify c, instead returns the mutated c.  Probably OK:  3/8/21
function mutate_chromosome!( c::Chromosome, funcs::Vector{Func}, mutate_location::Int64=0;
    insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  #println("mutate_chromosome! insert_gate_prob: ",insert_gate_prob)
  if rand() < insert_gate_prob && c.params.numinteriors < MyIntBits( MyInt )
    return(insert_gate!(c),false)
  elseif rand() < delete_gate_prob && c.params.numinteriors > c.params.numinputs 
    # Don't reduce number of gates less than number of interiors
    return (delete_gate!(c),true)  # active set to true, but has no meaning here
  end
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
    #println("func ",rand(1:100))
  elseif mutate_location <= num_funcs_to_mutate + sum(num_inputs_list)  # mutate an interior
    ml = mutate_location
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
    #print("ml: ",mutate_location,"  ")
    #print_circuit(c)
    #println("inputs: ",rand(1:100))
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
  set_active_to_false!(c)
  (c,active)
end

#= returns the component of the neutral set of the phenotype output_values(circuit) that contains circuit.
# This function is now defined in Neutral_component.jl
function neutral_component( circuit::Chromosome, funcs::Vector{Func}, component::Set=Set([]) )::Set
  ov = output_values(circuit)
  component = Set([])
  stack = Int64[]  # create a stack with Int64 elements
  chi = chromosome_to_int( circuit, funcs )
  push!( stack, chi )
  push!( component, chi )
  while length(stack) > 0
    cci = pop!(stack)
    cch = int_to_chromosome( cci, circuit.params, funcs )
    mcn = filter(cch->ov==output_values(cch), mutate_all( circuit, funcs, output_outputs=false, output_circuits=true ))
    mcn_pairs = map( c->(c,chromosome_to_int( c, funcs )), mcn )
    for ci_pair in mcn_pairs
      if !(ci_pair[2] in component)
        push!(component,ci_pair[2])
=#

# Mutates c in all possible ways. 
# The result depends on the keyword arguments.
# If output_outputs && !output_circuits (default) returns list of the outputs of chromsomes produced by mutation   
# If !output_outputs && output_circuits returns list of the chromsomes produced by mutation   
# If output_outputs && output_circuits returns both the list of outputs and the list of chromosomes as a pair
# If robustness_only==true returns the pair: (avg_robustness, evolvability)
# Deterministic if all functions in default_funcs() have the same arity
# TODO:  simplify by using the model of mutate_all() in LinChromosome.jl
function mutate_all( c::Chromosome, funcs::Vector{Func};
      robustness_only::Bool=false, output_outputs::Bool=true, output_circuits::Bool=false )
  #println("mutate_all: numlevelsback: ", c.params.numlevelsback )
  #sav_c = deepcopy(c)
  if robustness_only
    output_outputs = output_circuits = false
    robustness_sum = 0.0; robustness_count = 0
    result = Vector{MyInt}[]   # Save outputs to give a measure of evolvability
  elseif output_outputs && !output_circuits
    result = Vector{MyInt}[]
  elseif !output_outputs && output_circuits 
    result = Chromosome[]
  elseif output_outputs && output_circuits 
    #result = Tuple{Vector{MyInt},Chromosome}[]
    outputs_list = Vector{MyInt}[]
    ch_list = Chromosome[]
  end
  context = construct_context(c.params.numinputs)
  orig_output = output_values(c)
  num_mutate_locs = num_mutate_locations( c, funcs )
  #println("num_mutate_locs: ",num_mutate_locs)
  interiors_inputs_list = [ c.interiors[i].inputs for i = 1:c.params.numinteriors ]
  num_inputs_list = map(length, interiors_inputs_list)
  num_funcs_to_mutate = length(funcs) > 1 ? c.params.numinteriors : 0
  for mutate_location = 1:num_mutate_locs
    #println("mutate location: ",mutate_location)
    new_c = output_circuits ? deepcopy(c) : c
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
        new_c = deepcopy(c)  # This is necessary:  test done on 9/24/21
        #new_c = output_circuits ? deepcopy(c) : c
        new_c.interiors[mutate_location] = InteriorNode(new_func,new_inputs)
        deactivate_chromosome!(new_c)
        new_output = execute_chromosome(new_c,context)
        if robustness_only
          robustness_sum = (new_output == orig_output) ? robustness_sum + 1 : robustness_sum
          robustness_count+=1
          push!(result,new_output)
        elseif output_outputs && !output_circuits
          push!(result,new_output)
        elseif !output_outputs && output_circuits
          push!(result,new_c)
        elseif output_outputs && output_circuits
          #push!(result,(new_output,new_c))
          push!(outputs_list,new_output)
          push!(ch_list,new_c)
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
        new_c = deepcopy(c)   # Is this necessary???  Maybe only when the chromosome is returned
        new_c.interiors[i].inputs[j] = intindex
        #println("intindex: ",intindex,"  (i,j): ",(i,j),"  new_c.interiors[i].inputs: ",new_c.interiors[i].inputs)
        deactivate_chromosome!(new_c)
        new_output = execute_chromosome(new_c,context)
        if robustness_only
          robustness_sum = (new_output == orig_output) ? robustness_sum + 1 : robustness_sum
          robustness_count+=1
          push!(result,new_output)
        elseif output_outputs && !output_circuits
          push!(result,new_output)
        elseif !output_outputs && output_circuits
          push!(result,new_c)
        elseif output_outputs && output_circuits
          #push!(result,(new_output,new_c))
          push!(outputs_list,new_output)
          push!(ch_list,new_c)
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
    #println("length result: ",length(result))
  end
  if robustness_only
    (robustness_sum/robustness_count, length(unique(result))/length(result)) # pair of avg robustenss and evolvability
  elseif output_outputs && output_circuits
    return (outputs_list,ch_list)
  else
    result
  end
end

# Returns the list of chromosomes returned by mutate_all whose output is neutral, i. e., whose output is the phenotype of c
function mutate_all_neutral( c::Chromosome, funcs::Vector{Func} )
  ph = output_values( c )
  (outputs_list,ch_list) = mutate_all(c,funcs,output_circuits=true,output_outputs=true);
  ch_list[findall(x->x[1]==ph[1],outputs_list)];
end

# Supposed to return the number of ways to mutate c
# However, mutate_all returns more than this number of mutations.
function num_mutate_locations( c::Chromosome, funcs::Vector{Func} )
  if length(c.interiors) == 0
    return 0
  end
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
  #check_recursive(c)
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

# Counts number active nodes which includes inputs
#= Overwritten by another version below
function number_active( c::Chromosome )
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
  #if length(c.outputs)==0 || !c[c.outputs[1].input].active   # if chromosome has not been executed # Possible fix for rare bug
    context = construct_context(c.params.numinputs)
    execute_chromosome(c,context)
  end
  num_act_in = reduce(+,[c.inputs[i].active for i = 1:length(c.inputs)])
  num_act_int =  reduce(+,[c.interiors[i].active for i = 1:length(c.interiors)])
  num_act_in + num_act_int
end
=#

# Counts number active gates which does not include inputs
function number_active_gates( c::Chromosome )
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
    context = construct_context(c.params.numinputs)
    execute_chromosome(c,context)
  end
  num_act_int =  reduce(+,[c.interiors[i].active for i = 1:length(c.interiors)])
end

function number_active( c::Chromosome )
  number_active_gates( c )
end

# Requires that cache values are set in the various versions of evaluat_node()
function node_values( c::Chromosome )
  output_values = [c[c.outputs[i].input].cache for i = 1:length(c.outputs)]
  input_values = [c.inputs[i].cache for i = 1:length(c.inputs)]
  interior_values = [c.interiors[i].cache for i = 1:(length(c.interiors))]
  #output_values = [c.interiors[i].cache for i = (length(c.interiors)-length(c.outputs)+1):length(c.interiors)]
  (input_values,interior_values,output_values)
end

function output_values( c::Chromosome )
  if length( c.outputs ) == 0
    println( "empty chromosome in output_values()")
    return [ MyInt(0) ]
  end
  if !c[c.outputs[1].input].active   # if chromosome has not been executed
    context = construct_context(c.params.numinputs)
    return execute_chromosome(c,context)
  else
    return [c[c.outputs[i].input].cache for i = 1:length(c.outputs)]
  end
end

function output_values( c::Chromosome, funcs::Vector{Func} )
  output_values( c )
end

# Sets the active field of all input and interior nodes to be false
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

# Number of bits where x and y differ.  x and y are interpreted as bit strings
function hamming( x::MyInt, y::MyInt )
  xr = xor( x, y )
  result = 0
  my_one = MyInt(1)
  my_zero = MyInt(0)
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
function hamming_distance( x::Vector{MyInt}, y::Vector{MyInt}, c::Circuit )
  numinputs = c.params.numinputs
  @assert length( x ) == length( y )
  sum( hamming_distance( x[i], y[i], numinputs ) for i = 1:length(x))/length(x)
end

# Hamming distance between Goals which is between 0 and 1.
function hamming_distance( x::Vector{MyInt}, y::Vector{MyInt}, numinputs::Int64 )
  @assert length( x ) == length( y )
  sum( hamming_distance( x[i], y[i], numinputs ) for i = 1:length(x))/length(x)
end

# the minimum of the hamming distance and 1+hamming_distance( x, bitwise_complement(y))
function ihamming_distance( x::MyInt, y::MyInt, numinputs::Int64 )
  Ones = Main.CGP.construct_ones(numinputs)[numinputs]
  complement_y = xor(y,Ones)
  #@printf("0x0%0x\n",complement_y)
  #println("hd:  ",hamming_distance(x,y,numinputs))
  #println("chd: ",hamming_distance(x,complement_y,numinputs))
  #println("xhd: ",1.0/numinputs+hamming_distance(x,complement_y,numinputs))
  return min(hamming_distance(x,y,numinputs),1.0/2^numinputs+hamming_distance(x,complement_y,numinputs))
end

function ihamming_distance( x::Vector{MyInt}, y::Vector{MyInt}, numinputs::Int64 ) 
  sum( ihamming_distance( x[i], y[i], numinputs ) for i = 1:length(x))/length(x) 
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
# julia>circuit((1,2,3), ((4,OR,1,2), (5,AND,2,3), (6,XOR,4,5)),levsback=3)
#    creates a circuit with 3 input nodes and 3 gate nodes and 1 implicit output node.
#  The levsback keyword parameter is necessary to set the circuit parameters correctly.
# gate nodes have 4 fields: 
#    node_index:  should be the index of the node.  Not used in construcing the circuit, but checked for correctness
#    node_function:  See Func.jl for options
#    node_input1:  an integer in the range from 1 to the index of the node minus 1
#    node_input2:  an integer in the range from 1 to the index of the node minus 1
# Assumes nodearity==2 and numoutputs == 1
# Creates a Chromosome (circuit) from the concise format.
# Example:
# julia>circuit((1,2,3), ((4,OR,1,2), (5,AND,2,3), (6,XOR,4,5)),levsback=3)
#    creates a circuit with 3 input nodes and 3 gate nodes and 1 implicit output node.
#  The levsback keyword parameter is necessary to set the circuit parameters correctly.
# gate nodes have 4 fields:
#    node_index:  should be the index of the node.  Not used in construcing the circuit, but checked for correctness
#    node_function:  See Func.jl for options
#    node_input1:  an integer in the range from 1 to the index of the node minus 1
#    node_input2:  an integer in the range from 1 to the index of the node minus 1
# Example:
# julia>circuit((1,2,3), ((4,OR,1,2), (5,AND,2,3), (6,XOR,4,5)),levsback=3)
#    creates a circuit with 3 input nodes and 3 gate nodes and 1 implicit output node.
#  The levsback keyword parameter may be necessary to set the circuit parameters correctly.
# Example:
# julia>circuit((1,2,3), ((4,OR,1,2), (5,AND,2,3), (6,XOR,4,5)), (5,6), levsback=3)
#    creates a circuit with 3 input nodes and 3 gate nodes and output nodes 5 and 6
#  The levsback keyword parameter may be necessary to set the circuit parameters correctly.
# Assumes nodearity==2 and numoutputs == 1
function circuit( inputs::Tuple, gates::Tuple, outputs::Tuple=(); levsback::Int64=length(inputs)+length(gates))
  if outputs == ()
    outputs= (OutputNode( length(inputs)+length(gates )),)
    #println("outputs: ",outputs)
    #println("typeof outputs: ",typeof(outputs))
  else
    outputs = Tuple(OutputNode(i) for i in outputs)
  end
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
  p = Parameters( numinputs=length(inputs), numoutputs=length(outputs), numinteriors=length(gates), numlevelsback=levsback )
  in_nodes = [InputNode(i) for i in inputs]
  gate_nodes = [InteriorNode(g[2], [g[3],g[4]] ) for g in gates ]
  out_nodes = [OutputNode(i.input) for i in outputs ]
  c = Chromosome( p, in_nodes, gate_nodes, out_nodes, 0.0, 0.0 )
  if errors
    println("correct input to this function should have been: ")
    print_circuit(c)
  end
  c
end

#=
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
=#
  
# Outputs the concise format of Chromosome (circuit) c to IO stream f
# print_node_tuple() and gate_tuple() are defined in Node.jl
function print_circuit( f::IO, c::Chromosome; include_fitness::Bool=false, include_robustness::Bool=false,
    include_pheno::Bool=false ) 
  @assert c.params.nodearity==2
  @assert c.params.numinputs>1
  print(f,"circuit(")
  print_node_tuple(f, c.inputs )
  print(f,",")
  gate_tuple(f, c.interiors, c.params.numinputs )
  if c.params.numoutputs>1
    print(f,",")
    print_node_tuple(f, c.outputs )
  end
  if include_fitness
    @printf(f,", %5.3f",c.fitness)
  end
  if include_robustness
    if typeof(c.robustness) == PredType
      print(f,", ",c.robustness)
    elseif typeof(c.robustness) == Float64
      @printf(f,", %5.3f",c.robustness)
    end
  end
  if include_pheno   
    if c.params.numoutputs == 1
      @printf(f,", 0x%x",output_values(c)[1])
    else
      print(f,", [")
      for i = 1:(c.params.numoutputs-1)
        @printf(f,"0x%x, ",output_values(c)[i])
      end
      @printf(f,"0x%x]",output_values(c)[c.params.numoutputs])
    end
  end
  println(f,")")
end

function print_circuit( c::Chromosome; include_fitness::Bool=false, include_robustness::Bool=false, include_pheno::Bool=false ) 
  print_circuit( Base.stdout, c, include_fitness=include_fitness, include_robustness=include_robustness, include_pheno=include_pheno ) 
end    

function print_circuit( c::Chromosome, funcs::Vector{Func}; include_fitness::Bool=false, include_robustness::Bool=false, include_pheno::Bool=false ) 
  print_circuit( Base.stdout, c, include_fitness=include_fitness, include_robustness=include_robustness, include_pheno=include_pheno ) 
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


# Returns a Int128 integer which is unique to this chromsome.
# The values map one-to-one onto the integers in the range from 0 to count_circuits_ch(p,funcs)-1.
# Generates a Int128 integer corresponding to chromosome ch.  The integer is unique for parameters p.  
# Assumes that outputs come from the last numoutputs interior nodes and thus don't affect the number of chromosomes.
# See test/testChromosome.jl for a test----but only where all gates have arity maxarity.
function chromosome_to_int( ch::Chromosome, funcs::Vector{Func}=default_funcs(ch.params.numinputs); maxarity::Int64=2 )
  p = ch.params
  if count_circuits_ch( p, nfuncs=length(funcs) ) < Float64(typemax(Int128))
    result = Int128(0)
  else
    result = BigInt(0)
  end
  for i = 1:length(ch.interiors)
    (gint,multiplier) = gate_int( i+p.numinputs, ch, funcs )
    result = result*multiplier+gint
    #println("i: ",i,"  gate_int: ",gate_int( i+p.numinputs, ch, funcs ),"  multiplier: ",multiplier, "  result: ",result)
  end
  #println("ch_to_int result: ",result)
  result
end

# Same as chromosome_to_int() but with a different name
function circuit_to_int( ch::Chromosome, funcs::Vector{Func}=default_funcs(ch.params.numinputs); maxarity::Int64=2 )
  chromosome_to_int( ch, funcs, maxarity=maxarity )
end

# i is the Chromosome index of the gate as defined by getindex(ch,i) rather than the index in ch.interiors.
function gate_int( i::Int64, ch::Chromosome, funcs::Vector{Func} )
  interior_node = getindex(ch,i)
  #println("interior_node.func: ",interior_node.func)
  p = ch.params
  func_list = [ f.func for f in funcs ]   # necessary because the == operator doesn't work on structs
  func_int = findfirst(x->x==interior_node.func.func,func_list)[1]-1
  gate_inputs = interior_node.inputs
  all_possible_inputs = inputs_list( p.nodearity, max(1,i-p.numlevelsback), i-1 )
  #println("i: ",i,"  gate_inputs: ",gate_inputs,"  all_possible_inputs: ",all_possible_inputs)
  inputs_int = findfirst(x->x==gate_inputs,all_possible_inputs)-1
  multiplier = length(all_possible_inputs)*length(func_list) 
  (inputs_int*length(func_list)+func_int, multiplier)
end

# returns a list of the possilbe inputs to a chromosome where minval and maxval are 
#   determined by the parameters and the position of the gate in the chromosome.
# In all non-recursive calls to this function, numgateinputs is maxarity for the gates that we use in the evolvability paper.
# In recursive calls, numinputs is decreased
# Example for Parameters(2,1,4,3):  
#     inputs_list(2,2,4): [[2, 2], [2, 3], [2, 4], [3, 2], [3, 3], [3, 4], [4, 2], [4, 3], [4, 4]]
function inputs_list( numgateinputs::Int64, minval::Int64, maxval::Int64 )
  #println("inputs_list numgateinputs: ",numgateinputs,"  minval: ",minval,"  maxval: ",maxval )
  if numgateinputs == 1
    return [ [i] for i = minval:maxval ]
  end
  result = inputs_list( numgateinputs-1, minval, maxval )
  new_result = Vector{Int64}[]  # empty list of phenotypes1
  for r in result
    for i = minval:maxval
      dcr = deepcopy(r)
      push!(dcr,i)
      push!(new_result,dcr)
      #println("i: ",i,"  dcr: ",dcr)
    end
  end
  #println("inputs_list: numgateinputs: ",numgateinputs,"  new_result: ",new_result)
  new_result
end     

# Produces a list of all circuits corresponding the paramters p and funcs
function enumerate_circuits_ch( p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs); maxarity::Int64=2 )
  enumerate_circuits_ch( p, p.numinteriors, funcs, maxarity=maxarity )
end

# Recursive helper function
function enumerate_circuits_ch( p::Parameters, numints::Int64, funcs::Vector{Func}; maxarity::Int64=2 )
  #println("ec: numints: ",numints)
  if numints == 0
    result = [Chromosome( p, map(i->InputNode(i),collect(1:p.numinputs)), InteriorNode[], [OutputNode(p.numinputs)], 0.0, 0.0 )]
    return result
  end
  prev_result = enumerate_circuits_ch( p, numints-1, funcs ) 
  println("typeof(prev_result): ",typeof(prev_result),"  length(prev_result): ",length(prev_result))
  #println("prev_result: ",map(x->print_circuit(x),prev_result)
  result = Chromosome[]
  for prev_ch in prev_result
    #println("numints: ",numints,"  length(prev_result): ",length(prev_result))
    for inputs in inputs_list( maxarity, max(1,p.numinputs+numints-p.numlevelsback), p.numinputs+numints-1 )
      for func in funcs
        new_interiors = vcat(prev_ch.interiors,[InteriorNode(func,inputs)])
        new_ch = Chromosome( p, prev_ch.inputs,new_interiors,[OutputNode(p.numinputs+numints)],0.0,0.0) 
        push!(result,new_ch)
      end
    end
  end
  result
end

function int_to_chromosome( ch_int::Integer, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs); maxarity::Int64=2 )
  if count_circuits_ch( p, nfuncs=length(funcs) ) < Float64(typemax(Int128))
    ch_int = Int128(ch_int)
  else
    result = BigInt(ch_int)
  end
  nfuncs = length(funcs)
  inputnodes = [ InputNode(i) for i = 1:p.numinputs ]
  interiors = InteriorNode[]
  for i in p.numinteriors:-1:1
    inputslist = inputs_list( p.nodearity, max(1,p.numinputs+i-p.numlevelsback), i+p.numinputs-1 )
    fct = funcs[ch_int % nfuncs + 1 ]
    ch_int = div( ch_int, nfuncs )
    inputs = inputslist[ ch_int % length(inputslist) + 1 ]
    ch_int = div( ch_int, length(inputslist) )
    push!( interiors, InteriorNode( fct, inputs ) )
  end
  outputs = [ OutputNode(i) for i = (p.numinputs + p.numinteriors -p.numoutputs + 1):p.numinputs + p.numinteriors]
  Chromosome( p, inputnodes, reverse(interiors), outputs, 0.0, 0.0 )
end

# Same as int_to_chromosome() but with a different name
function int_to_circuit( ch_int::Integer, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs); maxarity::Int64=2 )
  int_to_chromosome( ch_int, p, funcs, maxarity=maxarity )
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

# Returns an equivalent chromosome where the inputs to every gate are in increasing order
function normalize_chromosome( ch::Chromosome )
  normalize_chromosome!(deepcopy(ch))
end

# Modifies chromosome ch so that the inputs to every gate are in increasing order
function normalize_chromosome!( ch::Chromosome )
  for i in ch.interiors
    if i.inputs[1] > i.inputs[2]
      tmp = i.inputs[1]; i.inputs[1] = i.inputs[2]; i.inputs[2] = tmp
    end
  end
  cah
end

# Returns the number of genotypes mapping to each phenotype for the parameters/funcs setting.
# This is for Chromosomes.
function count_genotypes_ch( p::Parameters, funcs::Vector{Func}=default_funcs(p) )::DataFrame
  num_circuits = Int(count_circuits_ch( p ))
  genotype_counts = zeros( Int64, 2^2^p.numinputs )
  for i = 1:num_circuits
    ov = output_values( int_to_chromosome(i,p,funcs) )[1]
    if ov+1 > 2^2^p.numinputs
      println("i: ",i,"  ov: ",ov)
    end
    genotype_counts[ov+1] += 1
  end
  #genotype_counts
  count_genotypes_table( p, funcs, genotype_counts )
end

# Returns the number of genotypes mapping to each phenotype for the parameters/funcs setting.
# This is for Chromosomes.
function count_genotypes_ch_mt( p::Parameters, funcs::Vector{Func}=default_funcs(p) )::DataFrame
  genotype_counts = [ Threads.Atomic{Int64}(0) for i= 1:2^2^p.numinputs ]
  num_circuits = Int(count_circuits_ch( p ))
  Threads.@threads for i = 1:num_circuits
    ov = output_values( int_to_chromosome(i,p,funcs))[1]
    Threads.atomic_add!( genotype_counts[ov+1], 1 )
  end
  gc = map( ph->ph[], genotype_counts )
  count_genotypes_table( p, funcs, gc )
end

# Return the number of circuits for parameters p and number of funcs nfuncs 
function count_circuits_ch( p::Parameters, funcs::Vector{Func} )::BigFloat
  gc = count_circuits_ch( p, nfuncs=length(funcs) )
end

# Return the number of chromosomes for parameters p and number of funcs nfuncs if nfuncs>0
# If nfuncs==0, nfuncs is reset to be length(default_funcs(p.numinputs))
function count_circuits_ch( p::Parameters; nfuncs::Int64=5 )::BigFloat
  @assert p.numoutputs == 1   # Not tested for more than 1 output, but probably works in this case.
  nfuncs = nfuncs==0 ? length(default_funcs(p.numinputs)) : nfuncs
  multiplier = BigFloat(1)
  k = 0
  for i = 1:p.numinteriors
    # uncomment for debugging
    #println("i: ",i,"  multiplier: ",multiplier,"  k: ",k)
    multiplier *= nfuncs
    for j = 1:p.nodearity
      k = min(p.numlevelsback,i-1+p.numinputs)
      multiplier *= k
      #println("i: ",i,"  j: ",j,"  multiplier: ",multiplier,"  k: ",k)
    end
    #println("i: ",i,"  nfuncs: ",nfuncs,"  k: ",k,"  multiplier: ",multiplier)
    #println("i: ",i,"  nfuncs: ",nfuncs,"  k: ",k,"  log multiplier: ",log10(multiplier))
  end
  multiplier^p.numoutputs
end

# 
function count_genotypes_table( p::Parameters, funcs::Vector{Func}, genotype_counts::Vector{Int64} )::DataFrame
  @assert length(genotype_counts) == 2^2^p.numinputs
  df = DataFrame(
    :phenotype=>map( ph->[MyInt(ph)], 0:(2^2^p.numinputs-1)),
    :count=>genotype_counts
  )
end

# Return the number of normalized circuits for parameters p and funcs 
function count_circuits_ch_normalize( p::Parameters, funcs::Vector{Func} )
  count_circuits_ch_normalize( p, nfuncs=length(funcs) )
end

# Return the number of normalized circuits for parameters p and number of funcs nfuncs 
function count_circuits_ch_normalize( p::Parameters, nfuncs::Int64 )
  count_circuits_ch_normalize( p, nfuncs=nfuncs )
end

# Return the number of normalized circuits for parameters p and number of funcs nfuncs if nfuncs>0
# If nfuncs==0, nfuncs is reset to be length(default_funcs(p.numinputs))
function count_circuits_ch_normalize( p::Parameters; nfuncs::Int64=0 )
  @assert p.numoutputs == 1   # Not tested for more than 1 output, but probably works in this case.
  nfuncs = nfuncs==0 ? length(default_funcs(p.numinputs)) : nfuncs
  multiplier = 1.0
  for i = 1:p.numinteriors
    n_inputs = min(p.numlevelsback,i-1+p.numinputs)  # Number of possible values for inputs to this interior node
    mult = (n_inputs+1)*n_inputs/2
    multiplier *= mult*nfuncs
    println("i: ",i,"  n_inputs: ",n_inputs,"  mult: ",mult,"  multiplier: ",multiplier)
  end
  multiplier
end

# Mutates chromosome c by inserting a gate or deleting a gate
function mutate_num_gates!( c::Chromosome; prob_insert::Float64=0.5 )
  if rand() <= prob_insert 
    c = insert_gate!(c) 
  else 
    c = insert_gate!(c) 
  end
  c
end

# Insert a gate into c so that output_values(c) is unchanged.
# The new gate is inserted at position new_gate_index.
# The following gates are shifted right by 1.
# The new gate replaces an input connection to the existing gate a position new_gate_index (before the shift)
# The output of the new gate is the input to the existing gate.
# One input to the new gate is the source of the input to the exisiting gate. 
# The other is chosen randomly according to the levelsback constraint.
function insert_gate!( c::Chromosome )
  funcs = default_funcs(c.params.numinputs)
  p = c.params
  # Chose a random interior node as the "new gate".
  # Don't insert a gate to replace an output gate:
  new_gate_index = rand(1:(p.numinteriors-p.numoutputs+1))  # Position of the new gate
  #println("new_gate_index: ",new_gate_index)
  p = c.params = Parameters( p.numinputs, p.numoutputs, p.numinteriors+1, p.numlevelsback+1 )
  new_gate = deepcopy(c.interiors[new_gate_index])  # will be modified below
  input_index = rand(1:2)
  new_gate.func = input_index == 1 ? IN1 : IN2  # IN1 and IN2 are gates that return input1 and input 2 respectively.
  c.interiors = vcat(c.interiors[1:(new_gate_index-1)],[new_gate],c.interiors[new_gate_index:end])
  for i = 1:p.numoutputs
    index = p.numinputs + p.numinteriors + i - p.numoutputs # use the last numoutputs interiors 
    c.outputs[i] = OutputNode(index) 
  end
  for i = (new_gate_index+1):p.numinteriors
    for j = 1:p.nodearity
      if c.interiors[i].inputs[j] >= p.numinputs + new_gate_index 
        c.interiors[i].inputs[j] += 1
      end
    end
  end
  # Modify one of the inputs of interiors[new_gate_index+1] to point to the new gate
  c.interiors[new_gate_index+1].inputs[input_index] = p.numinputs + new_gate_index
  set_active_to_false!(c)
  #println("insert_gate!  new interiors: ",c.params.numinteriors)
  c
end

function delete_gate!( c::Chromosome, interior_to_delete::Int64 )
  #println("delete_gate! interior_to_delete: ",interior_to_delete,"  active: ",c.interiors[interior_to_delete].active)
  p = c.params
  @assert interior_to_delete <= p.numinteriors - p.numoutputs
  gate_to_delete = c.interiors[interior_to_delete]
  new_levelsback = (p.numlevelsback > p.numinputs) ? p.numlevelsback-1 : p.numlevelsback  # Don't set numlevelsback to less than p.numinputs.
  p = c.params = Parameters( p.numinputs, p.numoutputs, p.numinteriors-1, new_levelsback )
  c.interiors = vcat(c.interiors[1:(interior_to_delete-1)],c.interiors[(interior_to_delete+1):end])
  for i = 1:p.numoutputs
    index = p.numinputs + p.numinteriors + i - p.numoutputs # use the last numoutputs interiors 
    c.outputs[i] = OutputNode(index) 
  end
  for i = interior_to_delete:p.numinteriors
    for j = 1:p.nodearity
      if c.interiors[i].inputs[j] >= interior_to_delete + p.numinputs
        c.interiors[i].inputs[j] -= 1
      end
    end
  end
  #println("delete_gate!  new interiors: ",c.params.numinteriors)
  c
end

function delete_gate!( c::Chromosome )
  output_values(c)  # needed to set active gates to active  
  if number_active_gates(c) == c.params.numinteriors
    dg = rand(1:(c.params.numinteriors-c.params.numoutputs))
    return delete_gate!( c, dg )
  end
  inactive_interior_list = filter!(i->!c.interiors[i].active, collect(1:c.params.numinteriors))
  dg = rand(inactive_interior_list)
  return delete_gate!( c, dg )
end

# The genotype distance between two chromosomes where the distance is the number of mismatches.
function geno_distance( ch1::Chromosome, ch2::Chromosome, funcs::Vector{Func} )
  p1 = ch1.params
  p2 = ch2.params
  @assert p1 == p2
  #gtuple1 = gate_tuple( p1, funcs, ch1.interiors )
  #gtuple2 = gate_tuple( p2, funcs, ch2.interiors )
  #mismatch = [ gtuple1[i] != gtuple2[i] ? 1 : 0 for i = 1:length(gtuple1) ]
  #glist1 = gate_list( p1, funcs, ch1.interiors )
  #glist2 = gate_list( p2, funcs, ch2.interiors )
  glist1 = gate_list( funcs, ch1.interiors )
  glist2 = gate_list( funcs, ch2.interiors )
  mismatch = [ glist1[i] != glist2[i] ? 1 : 0 for i = 1:length(glist1) ]
  reduce( +, mismatch )
end

# References the function combine_chromosomes( c1::Chromosome, c2::Chromosome ) in Chromosome.jl.
# Given a pair of phenotypes ph1 and ph2, test whether combining chromosomes c1 and c2 evolved separately to output ph1 and ph2 has
#   greater Tononi complexity than a chromosome c0 evolved to output the phenotype [ph1[1],ph2[1]].
function test_combine_complexity( p::Parameters, funcs::Vector{Func}, ph_pairs_list::Vector{Tuple{Goal,Goal}}, nreps::Int64, max_tries::Int64, max_steps::Int64; csvfile::String="" )
  df = DataFrame( :ph1=>Goal[], :ph2=>Goal[], :tcmplx_c1=>Float64[], :tkcmplx_c1=>Float64[], :kcmplx_c1=>Int64[], :tcmplx_c2=>Float64[], :tkcmplx_c2=>Float64[], :kcmplx_c2=>Int64[], 
      :cmplx_cmb=>Float64[], :cmplx_cc=>Float64[], :kcmplx=>Int64[], :tkcmplx=>Float64[] )
  p0 = Parameters( p.numinputs, 2*p.numoutputs, 2*p.numinteriors, p.numlevelsback )
  print_parameters(p0)
  for ph_pair in ph_pairs_list
    ph1 = ph_pair[1]
    result = kolmogorov_complexity( p, ph1, max_tries, max_steps, use_mut_evolve=false )
    kcmplx1 = result[2]  # Kolmogorov complexity
    tcmplx1 = result[4]  # Tononi complexity of minimal gate circuit
    ph2 = ph_pair[2]
    result = kolmogorov_complexity( p, ph2, max_tries, max_steps, use_mut_evolve=false )
    kcmplx2 = result[2]  # Kolmogorov complexity
    tcmplx2 = result[4]  # Tononi complexity of minimal gate circuit
    result = kolmogorov_complexity( p0, [ph1[1],ph2[1]], max_tries, max_steps, use_mut_evolve=false )
    kcmplx0 = result[2]  # Kolmogorov complexity
    tcmplx0 = result[4]  # Tononi complexity of minimal gate circuit
    for i = 1:nreps
      (cc,steps) = pheno_evolve( p0, funcs, [ph1[1],ph2[1]], 2*max_tries, max_steps )
      (c1,steps) = pheno_evolve( p, funcs, ph1, 2*max_tries, max_steps ) 
      (c2,steps) = pheno_evolve( p, funcs, ph2, 2*max_tries, max_steps ) 
      cmb = combine_chromosomes( c1, c2 )
      push!(df,(ph1,ph2,complexity5(c1),tcmplx1,kcmplx1,complexity5(c2),tcmplx2,kcmplx2,complexity5(cmb),complexity5(cc),kcmplx0,tcmplx0))
    end
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# Results test_combine_complxity() ")
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ",funcs)
      println(f,"# parameters p: ")
      print_parameters(f,p,comment=true)
      println(f,"# combined parameters p0: ")
      print_parameters(f,p0,comment=true)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end                  
  df
end

# Combines two chromosomes c1 and c2 into a chromosome whose outputs are the combined outputs of c2 and c1
function combine_chromosomes( c1::Chromosome, c2::Chromosome )
  @assert c1.params.numinputs == c2.params.numinputs
  p = Parameters( c1.params.numinputs, c1.params.numoutputs+c2.params.numoutputs, c1.params.numinteriors+c2.params.numinteriors, c1.params.numinteriors+c2.params.numinteriors+1 )
  inputs = deepcopy(c1.inputs)
  interiors1 = deepcopy(c1.interiors)
  interiors2 = deepcopy(c2.interiors)
  new_interiors = InteriorNode[]
  for i = 1:length(c2.interiors)
    new_interior1 = deepcopy(c1.interiors[i])
    new_interior2 = deepcopy(c2.interiors[i])
    new_interior1.active = false
    new_interior2.active = false
    for j = 1:length(new_interior.inputs)
      inc = new_interior.inputs[j] > p.numinputs ? c1.params.numinteriors : 0  # increment inputs that don't point to input nodes
      new_interior.inputs[j] += inc
    end
    push!(interiors,new_interior)
  end
  new_output2 = deepcopy(c2.outputs)
  for newout in new_output2
    newout.input = newout.input+c1.params.numinteriors
  end
  new_outputs = vcat(new_output2,deepcopy(c1.outputs))
  Chromosome( p, inputs, interiors, new_outputs, 0.0, 0.0 )
end

# Append the gates of chromosome c2 to those of chtomosome c1 and adjust parameters appropriately
# Gates of the second chromosone that refer to previous gates in that chromosome are adjusted appropriately
function append_chromosome( c1::Chromosome, c2::Chromosome )
  @assert c1.params.numinputs == c2.params.numinputs
  c2 = deepcopy(c2)
  p = Parameters( c1.params.numinputs, c2.params.numoutputs, c1.params.numinteriors+c2.params.numinteriors, c1.params.numinteriors+c2.params.numinteriors+1 )
  inputs = deepcopy(c1.inputs)
  set_active_to_false!( c1 )
  set_active_to_false!( c2 )
  new_interiors = deepcopy(c1.interiors)
  for i = 1:length(c2.interiors)
    println("i: ",i,"  c2.interiors[i]: ",c2.interiors[i])
    for j = 1:p.nodearity
      if c2.interiors[i].inputs[j] > c2.params.numinputs
        c2.interiors[i].inputs[j] += c1.params.numinteriors
      end
    end
    push!(new_interiors,deepcopy(c2.interiors[i]))
  end
  new_outputs = deepcopy(c2.outputs)
  for newout in new_outputs
    newout.input = newout.input+c1.params.numinteriors
  end
  Chromosome( p, inputs, new_interiors, new_outputs, 0.0, 0.0 )
end

function interleave_chromosomes( c1::Chromosome, c2::Chromosome )
  function increment_input( i::Int64, odd::Bool )
    if i > c1.params.numinputs
      #2*(i-c1.params.numinputs) + (odd ? 0 : 1) + c1.params.numinputs
      2*i - c1.params.numinputs + (odd ? 0 : 1) - 1
    else
      i
    end
  end
  @assert length(c1.interiors) == length(c2.interiors)
  new_interiors = InteriorNode[]
  for i = 1:length(c1.interiors)
    new_inputs = map( i->increment_input(i,true), c1.interiors[i].inputs )  
    new_interior = deepcopy(c1.interiors[i])  
    new_interior.inputs = new_inputs
    new_interior.active = false
    new_interior.cache = MyInt(0)
    push!( new_interiors, new_interior )
    new_inputs = map( i->increment_input(i, false), c2.interiors[i].inputs )  
    new_interior = deepcopy(c2.interiors[i])  
    new_interior.inputs = new_inputs
    new_interior.active = false
    new_interior.cache = MyInt(0)
    push!( new_interiors, new_interior )
  end
  cnew = deepcopy(c1)
  cnew.interiors = new_interiors
  ninputs = c1.params.numinputs
  cnew.params = Parameters( c1.params.numinputs, 2, c1.params.numinteriors+c2.params.numinteriors, c1.params.numlevelsback+c2.params.numlevelsback )  # increase numoutputs to 2
  cnew.outputs = OutputNode[ OutputNode(ninputs+2*length(c1.interiors)-1), OutputNode(ninputs+2*length(c1.interiors)) ]
  cnew
end

# Creates a chromosome whose interiors are concatenation of the interiors of c1 and c2
function chromosome_cat( c1::Chromosome, c2::Chromosome )
  @assert c1.params.numinputs == c2.params.numinputs
  @assert c1.params.numoutputs == c2.params.numoutputs == 1
  cnew = deepcopy(c1)
  cnew.params = Parameters( c1.params.numinputs, c1.params.numoutputs, 
      c1.params.numinteriors+c2.params.numinteriors, c1.params.numlevelsback+c2.params.numlevelsback )
  cnew.interiors = append!(cnew.interiors,c2.interiors)
  cnew.outputs = OutputNode[ OutputNode( cnew.params.numinputs + cnew.params.numinteriors ) ]
  cnew
end 

# Creates a new chromosome which is the same as ch except that it has another input---which does not contribute to the output
# However, the output changes mostly because the context changes.
function chromosome_add_input( ch::Chromosome )
  cnew = deepcopy( ch )
  set_active_to_false!(cnew)
  cnew.params = Parameters( ch.params.numinputs+1, ch.params.numoutputs, ch.params.numinteriors, ch.params.numlevelsback )
  cnew.inputs = push!( cnew.inputs, InputNode( cnew.params.numinputs ) )
  # adjust interior nodes to not refer to the new input node
  for intnode in cnew.interiors
    for i = 1:length(intnode.inputs)
      if intnode.inputs[i] > length(cnew.inputs)
        intnode.inputs[i] += 1
      end
    end
  end
  cnew.outputs = [ OutputNode( cnew.outputs[i].input+1 ) for i = 1:length(cnew.outputs) ]
  cnew
end

# adds a multiplexer at the end of an interleaved chromosome
function chromosome_add_multiplexer( ch::Chromosome )
  @assert ch.params.numinputs == 2
  @assert ch.params.numoutputs == 2
  x = ch.params.numinteriors
  println("x: ",x)
  output_gates = map( x->x.input, ch.outputs )
  cnew = chromosome_add_input( ch )
  mux = circuit((1,2,3), ((4,XOR,1,2), (5,AND,4,3), (6,XOR,2,5)))
  push!(cnew.interiors,InteriorNode(XOR,[output_gates[1],output_gates[2]]))
  push!(cnew.interiors,InteriorNode(AND,[3,length(cnew.interiors)+2]))
  push!(cnew.interiors,InteriorNode(XOR,[output_gates[2],length(cnew.interiors)+2]))
  cnew.outputs=OutputNode[OutputNode(length(cnew.interiors)+2)] 
  cnew.params = Parameters( ch.params.numinputs, 1, ch.params.numinteriors+3, ch.params.numlevelsback+3 )
  cnew
end
#=
function chromosome_add_multiplexer( ch::Chromosome )
  @assert ch.params.numinputs == 2
  x = length(ch.params.numinteriors)
  c1 = chromosome_add_input( ch )
  mux = circuit((1,2,3), ((4,XOR,1,2), (5,AND,4,3), (6,XOR,2,5)))
  push!(c1.interiors,InteriorNode(XOR,[x-1,x]))
  push!(c1.interiors,InteriorNode(AND,[3,x+1]))
  push!(c1.interiors,InteriorNode(XOR,[x,x+2]))
  c1
end
=#

# find the average relatioship between the number of one bits in a phenotype and its log redundancy, K complexity
# Written on 12/23/22.  Results were disappointing
function count_ones_properties( p::Parameters, funcs::Vector{Func} )
  @assert p.numinputs <= 4
  #println("methods robustness: ",methods(robustness) )
  kdict = kolmogorov_complexity_dict( p, funcs )
  rdict = redundancy_dict( p, funcs )
  pheno_list = collect(MyInt(0):MyInt(2^2^p.numinputs-1))
  count_ones_list = map( ph->count_ones(ph), pheno_list )
  lg_redund_list = map( ph->lg10(rdict[ph]), pheno_list )
  k_complexity_list = map( ph->kdict[ph], pheno_list )
  #robustness_list = map( ph->robustness( ph, funcs ), pheno_list )
  #DataFrame( :pheno=>pheno_list, :log_redund=>lg_redund_list, :k_complexity=>k_complexity_list, :robustness=>robustness_list )
  DataFrame( :pheno=>pheno_list, :count_ones=>count_ones_list, :log_redund=>lg_redund_list, :k_complexity=>k_complexity_list )
end

#= Not correct or finished
function redundancy_density( p::Parameters, funcs::Vector{Func}, nsamples::Int64 )
  rdict = redundancy_dict( p, funcs )
  kdict = kolmogorov_complexity_dict( p, funcs )
  gdf = read_dataframe( "../data/counts/count_outputs_ch_5funcs_4inputs_10gates_5lb_EG.csv" )
  glogredund = map( r->lg10(r), gdf.ints10_5 )
end
=#

