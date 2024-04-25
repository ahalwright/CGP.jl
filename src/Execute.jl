export execute_chromosome, concatenate_outputs, create_count_function_list, increment_count_function_list, print_function_list
export evaluate_node, create_count_function_hash, increment_count_function_hash, print_count_function_hash, print_count_function_summary
export execute_chromosome_ft, evaluate_node_ft, execute_chromosome_fwd, eval_interior
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
    node.active = true
    #println("eval input node: ",node)
    #print("context[node.index]:")
    #Printf.@printf("0x0%4x\n",context[node.index]) 
    #=
    if node.active
      return node.cache
    end
    node.cache = context[node.index]
    =#
    return context[node.index]
end

# Restored caching on 10/7/21
function evaluate_node(c::Chromosome, node::InteriorNode, context::Vector)
    #println("en eval interior node: ",node,"  contxt: ",context)
    prev_cache = node.cache
    func = node.func
    if node.active
      return node.cache
    end
    #println("node.inputs: ",node.inputs)
    args = map(node.inputs[1:func.arity]) do index
        #println("index: ",index)
        evaluate_node(c, c[index], context)
    end
    #args = map(x->context[x],node.inputs)
    #println("args: ",args)
    result = CGP.apply(func.func, args) 
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
# context should be the context for numinputs.
# Not currently used.  4/22/24
function eval_interior( node::InteriorNode, context::Vector )
    println("ei eval interior node: ",node)
    func = node.func
    args = [context[i] for i = 1:func.arity]
    val = Main.CGP.apply(func.func, args)
    println("args: ",args,"  val: ",val)
    val
end

function evaluate_node(c::Chromosome, node::OutputNode, context::Vector)
    index = node.input
    result = evaluate_node(c, c[index], context)
    #println("eval output node: ",node,"  index: ",index,"  result: ",result)
    #result = evaluate_node(c, c[index], context)
    return result
end

# Execute chromosme c using context for the inputs and return outputs
# c is modified in that the active and cache fields are set
# This is the function used by the output_values() function
function execute_chromosome(c::Chromosome, context::Vector; permute_context::Bool=false )
    #println("Executing chromosome")
    #print_chromosome(c)
    #println("context: ",context)
    return [evaluate_node(c, node, context) for node = c.outputs]
end

# ft means fault_tolerance.  Not currently used, might be deleted
function evaluate_node_ft(c::Chromosome, node::OutputNode, context::Vector, perturb_index::Int64 )
    index = node.input
    #println("ft eval output node: ",node,"  index: ",index)
    if index == perturb_index
        #println("ft out index == perturb_index")
        result = NOT.func(evaluate_node_ft(c, c[index], context,perturb_index))
    else
       result = evaluate_node_ft(c, c[index], context,perturb_index)
    end
    #print("ft out index: ",index)
    #Printf.@printf("ft out result: 0x%2x\n",result)
    return result
end

function evaluate_node_ft(c::Chromosome, node::InputNode, context::Vector, perturb_index::Int64 )
    #println("ft eval input node: ",node)
    #Printf.@printf("ft result:  0x0%4x\n",context[node.index]) 
    return context[node.index]
end

function evaluate_node_ft(c::Chromosome, node::InteriorNode, context::Vector, perturb_index::Int64 )
    #println("ft eval interior node: ",node)
    func = node.func
    args = map(node.inputs[1:func.arity]) do index
        #println("index: ",index)
        if index == perturb_index
          #println("  ft int index == perturb_index")
          NOT.func(evaluate_node_ft(c, c[index], context,perturb_index))
        else
          evaluate_node_ft(c, c[index], context,perturb_index)
        end
    end
    result = apply(func.func, args) 
    #print("ft int func: ",func,"  args: ",args )
    #Printf.@printf("  ft int result: 0x%2x\n",result)
    return result
end

# Fault tolerant version where the output of node c[perturb_index] is negated
# Note: ignores cache values for nodes
function execute_chromosome_ft(c::Chromosome, context::Vector, perturb_index::Int64 )
    #println("FT Executing chromosome")
    #print_chromosome(c)
    return [evaluate_node_ft(c, node, context,perturb_index) for node = c.outputs]
end

# Executes chromosome c by evaluating interior nodes in the order of the chromosome
# The default is to recursively execute notes starting with the output node(s)
function execute_chromosome_fwd(c::Chromosome, context::Vector )
  values = MyInt[]
  for i = 1:(c.params.numinputs + c.params.numinteriors + c.params.numoutputs)
    if i <= c.params.numinputs
      push!(values,context[i])
    elseif i <= c.params.numinputs + c.params.numinteriors 
      args = map(j->values[j],c[i].inputs)
      push!(values, eval_interior( c[i], args ) )
      #println("i: ",i,"  values: ",values)
    end
  end
  values[end]
end
#=
function print_function_list( function_list::Vector{MyFunc} )
  for i = convert(MyFunc,1):length(function_list)
    if function_list[i] > 0
      # The correct way to print when MyInt=UInt16, but wont work when MyInt is something else
      #Printf.@printf("0x%04x  %d\n",i-1,function_list[i])  
    end
  end
end
=#

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

function concatenate_outputs( numinputs::Integer, outputs::Vector{MyInt}  )
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
    #=
    if len <= length_counts_list_max_for_printing
      Printf.@printf(f,"0x%04x   %d\n",c[1],c[2])
    end
    =#
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

function poset( ch::Chromosome, funcs::Vector{Func} )
  p = ch.params
  len = length(ch.inputs) + length(ch.interiors) + length(ch.outputs)
  adj = zeros(Bool,len,len)
  for i = 1:len
    adj[i,i] = true
  end
  for i = len:-1:(len-length(ch.outputs)+1)
    adj[ch[i].input,i] = true
  end
  for i = (len-length(ch.outputs)):-1:(len-length(ch.outputs)-length(ch.interiors)+1)
    for j = 1:p.nodearity
      adj[ch[i].inputs[j],i] = true
    end
  end
  adj
end

# Warshall's algorithm
function transitive_closure!( adj::Matrix{Bool} )
  len = size(adj)[1]
  for i = 1:len
    for j = 1:len
      for k = 1:len
        if !adj[j,k]
          adj[j,k] = adj[j,i] && adj[i,k]
        end
      end # for k
    end # for j
  end # for i
end

function topo_sort_all( adj::Matrix{Bool} )
  len = size(adj)[1]
  path_list = Vector{Int64}[]
  topo_sort_recurse_all( adj, path_list, collect(1:len) )
end

function topo_sort_recurse_all( adj::Matrix{Bool}, path_list::Vector{Vector{Int64}}, map::Vector{Int64} )
  len = size(adj)[1]
  println("topo_sort_recurse: len: ",len)
  if len == 0
    return path
  end
  indegree = [ sum(adj[:,i])-1 for i = 1:len ]
  fmin_list = findall( x->x==0, indegree )
  if indegree[fmin_list[1]] == 0
    r = fmin_list[1]
    println("r: ",r,"  map: ",map)
    push!(path_list[1],map[r])
    deleteat!( map, r )
    println("path_list[1]: ",path_list[1],"  map: ",map)
    new_adj = adj[ 1:end .!= r, 1:end .!= r ]
    topo_sort_recursei_all( new_adj, path_list, map )
  end
  path_list
end

function topo_sort( adj::Matrix{Bool} )
  len = size(adj)[1]
  path = Int64[]
  topo_sort_recurse( adj, path, collect(1:len) )
end

function topo_sort_recurse( adj::Matrix{Bool}, path::Vector{Int64}, map::Vector{Int64} )
  len = size(adj)[1]
  println("topo_sort_recurse: len: ",len)
  if len == 0
    return path
  end
  indegree = [ sum(adj[:,i])-1 for i = 1:len ]
  fmin = findmin( indegree )
  if fmin[1] == 0
    r = fmin[2]
    println("r: ",r,"  map: ",map)
    push!(path,map[r])
    deleteat!( map, r )
    println("path: ",path,"  map: ",map)
    new_adj = adj[ 1:end .!= r, 1:end .!= r ]
    topo_sort_recurse( new_adj, path, map )
  end  
  path
end    

function topo_sort_check( adj::Matrix{Bool} )
  len = size(adj)[1]
  path = topo_sort( adj )
  chk_adj = adj[path,path]
  for i = 2:len
    for j = 1:(i-1)
      #println("(i,j): ",(i,j))
      if chk_adj[i,j] 
        return (i,j)
      end
    end
  end
  return "pass"
end
#
function remove_vertex( adj::Matrix{Bool}, r::Int64 )
  new_adj = adj[ 1:end .!= r, 1:end .!= r ]
end

function find_topo_sorts( adj::Matrix{Bool} )
  len = size(adj)[1]
  indegree = [ sum(adj[:,i])-1 for i = 1:len ]
  path = Int64[]
  discovered = falses( len )
  find_topo_sorts_recurse( adj, indegree, path, discovered )
end
    
function find_topo_sorts_recurse( adj::Matrix{Bool}, indegree::Vector{Int64}, path::Vector{Int64}, discoverd::Vector{Bool} )
  topo_sorts = Vector{Int64}[]
  for v = 1:len
    if indegree[v] == 0 && !discovered[v]
      for u = 1:len
        if adj[v,u] 
          indegree[u] -= 1
        end 
      end
      push!(path,v)
      discovered[v] = true
      find_topo_sorts_recurse( adj, indegree, path, discovered )
      for u = 1:len
        if adj[v,u] 
          indegree[u] += 1
        end 
      end
      pop!(path)
      discovered[v] = false
    end
  end
end
      
      
