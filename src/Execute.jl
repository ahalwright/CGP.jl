export execute_chromosome, concatenate_outputs, create_count_function_list, increment_count_function_list, print_function_list
export evaluate_node, create_count_function_hash, increment_count_function_hash, print_count_function_hash, print_count_function_summary
using Printf

# Limits Funcs to 4 inputs
# TODO:  Write this function generally
function apply( f::Function,args)
  if length(args) == 0
    f()
  elseif length(args) == 1
    f(args[1])
  elseif length(args) == 2
    f(args[1],args[2])
  elseif length(args) == 3
    f(args[1],args[2],args[3])
  elseif length(args) == 4
    f(args[1],args[2],args[3],args[4])
  end
end

function evaluate_node(c::Chromosome, node::InputNode, context::Vector)
    #println("eval input node: ",node)
    #print("context[node.index]:")
    #Printf.@printf("0x0%4x\n",context[node.index]) 
    if node.active
      return node.cache
    end
    node.active = true
    node.cache = context[node.index]
    return context[node.index]
end

function evaluate_node(c::Chromosome, node::InteriorNode, context::Vector)
    #println("eval interior node: ",node)
    prev_cache = node.cache
    func = node.func
    if node.active
      return node.cache
    end
    args = map(node.inputs[1:func.arity]) do position
        index = position
        #println("position: ",position)
        evaluate_node(c, c[index], context)
    end
    #println("args: ",args)
    result = apply(func.func, args) 
    #result = 0
    #print("func: ",func,"  args: ",args,"  result: ")
    #Printf.@printf("0x%2x\n",result)
    if node.active 
      @assert result == node.cache 
    end
    node.active = true
    node.cache = result
    return result
end

# Evaluate an interior node in isolation using context as the input
# context should be the context for nodearity.
# For debugging 
function eval_interior( node::InteriorNode, context::Vector )
    #println("eval interior node: ",node)
    func = node.func
    args = [context[i] for i = 1:func.arity]
    #println("args: ",args)
    Main.CGP.apply(func.func, args)
end

function evaluate_node(c::Chromosome, node::OutputNode, context::Vector)
    #println("eval output node: ",node)
    index = node.input
    return evaluate_node(c, c[index], context)
end

function execute_chromosome(c::Chromosome, context::Vector)
    #println("Executing chromosome")
    #print_chromosome(c)
    return [evaluate_node(c, node, context) for node = c.outputs]
end

function print_function_list( function_list::Vector{MyFunc} )
  for i = convert(MyFunc,1):length(function_list)
    if function_list[i] > 0
      Printf.@printf("0x%0x  %d\n",i-1,function_list[i])
    end
  end
end

# The purpose of the remaining functions in this file is to
#      generate random circuits and count the logic functions that are computed.
# Similar to results of Raman and Wagner 2011 to generate their Figure 2b.      
function create_count_function_list( numinputs::Integer, numoutputs::Integer )
  return fill( convert(MyFunc,0), 2^(numoutputs*2^numinputs ))
end

function increment_count_function_list( funct::MyFunc, function_list::Vector{MyFunc} )
  function_list[funct+1] += 1
end

function create_count_function_hash( numinputs::Integer, numoutputs::Integer )
  Dict{MyFunc,Int64}()
end

function increment_count_function_hash( funct::MyFunc, counts::Dict{MyFunc,Int64} )
  counts[funct] = get!( counts, funct, 0 ) + 1
end

function concatenate_outputs( numinputs::Integer, outputs::Vector{Main.CGP.MyInt}  )
  shift = 2^numinputs
  result = convert(MyFunc,outputs[1])
  for i = 2:length(outputs)
    result = convert(MyFunc,outputs[i]) | result << 2^numinputs
  end
  result
end
 
function print_count_function_hash( f::IOStream, counts::Dict{MyFunc,Int64} )
  length_counts_list_max_for_printing = 50
  println(f,"function count")
  counts_list = [(k,counts[k]) for k in keys(counts) ]
  sort!( counts_list, by=x->x[2], rev=true )
  len = length(counts_list)
  sum = 0
  for c in counts_list[1:len]
    if len <= length_counts_list_max_for_printing
      Printf.@printf(f,"0x%04x   %d\n",c[1],c[2])
    end
    sum += c[2]
  end
  println(f,"length counts list: ",len,"  sum: ",sum)
end

function print_count_function_hash( counts::Dict{MyFunc,Int64} )
  print_count_function_hash( Base.stdout, counts )
end

function print_count_function_summary( f::IOStream,  counts::Dict{MyFunc,Int64} )
  length_print_count_max_for_printing = 150
  println(f,"count function hash: ")
  counts_list = [(k,counts[k]) for k in keys(counts) ]
  sort!( counts_list, by=x->x[2] )
  current_count = counts_list[1][2]
  accum_count = 1
  accum_sum = 0
  print_count = 0
  println(f,"funct freq  count")
  for i = 2:length(counts_list)
    if counts_list[i][2] == current_count
      accum_count += 1
    else
      if print_count <= length_print_count_max_for_printing 
        println(f,"     ",current_count,":      ",accum_count)
      end
      print_count += 1
      accum_sum += current_count*accum_count
      current_count = counts_list[i][2]
      accum_count = 1
    end
  end
  accum_sum += current_count*accum_count
  if print_count <= length_print_count_max_for_printing 
    println(f,"     ",current_count,":      ",accum_count)
  end
  println(f,"print_count: ",print_count,"  accum_sum: ",accum_sum)
end

function print_count_function_summary( counts::Dict{MyFunc,Int64} )
  print_count_function_summary( f, counts )
end
