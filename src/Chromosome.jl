import Base.getindex

export Chromosome, print_chromosome, getindex, random_chromosome, mutate_chromosome, hamming

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
            c[inputs[i]].active = true
        end
        c.interiors[index - p.numinputs] = InteriorNode(func, inputs)
    end

    minindex = max(1,p.numinputs + p.numinteriors - p.numlevelsback + 1)
    maxindex = p.numinputs + p.numinteriors
    #println("output (minindex, maxindex): ",(minindex,maxindex))
    for i = 1:length(c.outputs)
        index = rand(minindex:maxindex)
        c.outputs[i] = OutputNode(index)
        c[index].active = true
    end
    return c
end

function mutate_chromosome( c::Chromosome, funcs::Vector{Func} )
  interiors_inputs_list = [ c.interiors[i].inputs for i = 1:c.params.numinteriors ]
  outputs_input_list = [ c.outputs[i].input for i = 1:c.params.numoutputs ]
  println("int_list: ", interiors_inputs_list, "out_list: ",outputs_input_list)
  num_inputs_list = map(length, interiors_inputs_list)
  num_inputs = sum(num_inputs_list) + c.params.numoutputs
  println("num_inputs_list: ",num_inputs_list,"  num_inputs: ",num_inputs)
  num_mutate_locations = c.params.numinteriors + num_inputs # Note: c.params.numinteriors = number funcs
  mutate_location = rand(1:num_mutate_locations)
  println("num_mutate_locations: ",num_mutate_locations,"  mutate_location: ",mutate_location)
  if mutate_location <= c.params.numinteriors  # mutate a func
    println("mutate a func")
    if length(funcs) > 1   # can't mutate func if there is only one possible func
      new_func_index = rand(1:length(funcs))
      println("new_funct_index: ",new_func_index,"  funcs[new_func_index]: ",funcs[new_func_index])
      println("  c.interiors[mutate_location].func: ",c.interiors[mutate_location].func)
      println("condition: ", funcs[new_func_index] == c.interiors[mutate_location].func )
      while( funcs[new_func_index] == c.interiors[mutate_location].func )
        new_func_index = rand(1:length(funcs))
        println("while new_funct_index: ",new_func_index)
      end
      println("new_func_index: ",new_func_index)
      c.interiors[mutate_location].func = funcs[new_func_index]
    end
  elseif mutate_location <= c.params.numinteriors + sum(num_inputs_list)  # mutate an interior
    println("mutate an interior")
    interior_mutate_location = mutate_location - c.params.numinteriors
    println("interior mutate_location: ",interior_mutate_location)
    i = 1  # i will be the index of the interior node one of whose inputs will be mutated
    num_inputs_sum = num_inputs_list[1]
    println("i: ",i,"  num_inputs_sum: ",num_inputs_sum)
    while(num_inputs_sum < interior_mutate_location)
      i += 1
      num_inputs_sum += num_inputs_list[i]
      println("i: ",i,"  num_inputs_sum: ",num_inputs_sum)
    end
    j  = num_inputs_sum - interior_mutate_location + 1  # j is the index of the input to mutate
    println("(i,j): ",(i,j),"  num_inputs_sum: ",num_inputs_sum)
    maxindex = i + c.params.numinputs -1
    minindex = max(1, i + c.params.numinputs - c.params.numlevelsback)
    intindex = rand(minindex:maxindex)
    if maxindex > minindex
      while( intindex == c.interiors[i].inputs[j] )
        intindex = rand(minindex:maxindex)
      end
    end 
    println("minindex: ",minindex,"  maxindex: ",maxindex,"  intindex: ",intindex)
    c.interiors[i].inputs[j] = intindex
    println("(i,j): ",(i,j),"  c.interiors[i].inputs: ",c.interiors[i].inputs)
  elseif mutate_location > c.params.numinteriors + sum(num_inputs_list)   # mutate an output
    println("mutate the output: ",mutate_location - c.params.numinteriors - sum(num_inputs_list ))
    minindex = max(1, c.params.numinteriors + c.params.numinputs - c.params.numlevelsback + 1 ) 
    maxindex = c.params.numinteriors + c.params.numinputs
    println("minindex: ",minindex,"  maxindex: ",maxindex)
    if maxindex > minindex   # no change if there is only one alternative for the output node index
      outindex = rand(minindex:maxindex)
      println("outindex: ",outindex)
      while outindex == c.outputs[mutate_location - c.params.numinteriors - sum(num_inputs_list)].input
        outindex = rand(minindex:maxindex)
      end
      c.outputs[mutate_location - c.params.numinteriors - sum(num_inputs_list)].input = outindex
    end
  end
  c
end

function print_chromosome( c::Chromosome )
  for i = 1:(c.params.numinputs + c.params.numinteriors + c.params.numoutputs)
    println(c[i])
  end
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

