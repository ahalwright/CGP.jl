import Base.getindex
export Chromosome, print_chromosome, getindex, random_chromosome, mutate_chromosome!, hamming
export num_mutate_locations, set_active_to_false, fraction_active, check_recursive, node_values
export number_active, number_active_old

mutable struct Chromosome
    params::Parameters
    inputs::Vector{InputNode}
    interiors::Vector{InteriorNode}
    outputs::Vector{OutputNode}
    fitness::Real
end

function Chromosome(p::Parameters)
    inputs = Array{InputNode}( undef, p.numinputs)
    interiors = Array{InteriorNode}( undef, p.numinteriors )
    outputs = Array{OutputNode}( undef, p.numoutputs)
    fitness = 0.0
    return Chromosome(p, inputs, interiors, outputs, fitness)
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
        inputs = Array{Integer}(undef, func.arity)
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
        println("output ",i," index: ",index)
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
function mutate_chromosome!( c::Chromosome, funcs::Vector{Func}, mutate_location::Integer=0 )
  num_mutate_locs = num_mutate_locations( c, funcs )
  if mutate_location == 0 
    mutate_location = rand(1:num_mutate_locs)
  end
  #println("num_mutate_locations: ",num_mutate_locs,"  mutate_location: ",mutate_location)
  interiors_inputs_list = [ c.interiors[i].inputs for i = 1:c.params.numinteriors ]
  num_inputs_list = map(length, interiors_inputs_list)
  num_funcs_to_mutate = length(funcs) > 1 ? c.params.numinteriors : 0
  if mutate_location <= num_funcs_to_mutate   # mutate a func
    #println("mutate a func ml:",mutate_location)
    #active = c[mutate_location].active
    new_func_index = rand(1:length(funcs))
    #println("new_funct_index: ",new_func_index,"  funcs[new_func_index]: ",funcs[new_func_index])
    #println("  c.interiors[mutate_location].func: ",c.interiors[mutate_location].func)
    #println("condition: ", funcs[new_func_index] == c.interiors[mutate_location].func )
    while( funcs[new_func_index] == c.interiors[mutate_location].func )
      new_func_index = rand(1:length(funcs))
      #println("while new_funct_index: ",new_func_index)
    end
    #println("new_func_index: ",new_func_index)
    old_func_arity = c.interiors[mutate_location].func.arity
    old_inputs = c.interiors[mutate_location].inputs
    new_func = funcs[new_func_index]
    new_inputs = Array{Integer}(undef, new_func.arity)
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
    #active = c[i].active
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
    println("mutate the output: ",mutate_location - num_funcs_to_mutate - sum(num_inputs_list ))
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
  c
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
  num_act_in = reduce(+,[c.inputs[i].active for i = 1:length(c.inputs)])
  num_act_int =  reduce(+,[c.interiors[i].active for i = 1:length(c.interiors)])
  num_act_in + num_act_int
end

function node_values( c::Chromosome )
  input_values = [c.inputs[i].cache for i = 1:length(c.inputs)]
  interior_values = [c.interiors[i].cache for i = 1:(length(c.interiors)-length(c.outputs))]
  output_values = [c.interiors[i].cache for i = (length(c.interiors)-length(c.outputs)+1):length(c.interiors)]
  (input_values,interior_values,output_values)
end


function fraction_active( c::Chromosome )
  number_active(c)/(c.params.numinputs+c.params.numinteriors)
end

function getindex(c::Chromosome, index::Integer)
    if index <= c.params.numinputs
        return c.inputs[index]
    end

    if index > c.params.numinteriors + c.params.numinputs
        return c.outputs[index - c.params.numinteriors - c.params.numinputs]
    end

    return c.interiors[index-c.params.numinputs]
end

function hamming( x::MyInt, y::MyInt )
  xr = xor( x, y )
  result = 0
  my_one = convert(MyInt,1)
  while xr != convert(MyInt,0)
    result +=  xr & my_one
    xr >>= 1
  end
  result
end

