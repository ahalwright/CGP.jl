
# returns the adjacency matrix corresponding to the partial order induce by genotype ch.
# adj[i,j] == true  if one of the inputs of node j is node i.
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
# Applying transitive_closure! to an adjacency matrix does not seem to affect the results of topo_sort_all()
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

# Finds one topological sort for adjacency matrix adj
# A topological sort is an total ordering (permutation) of the vertices compatible with the partial order given by adj.
function topo_sort( adj::Matrix{Bool} )
  len = size(adj)[1]
  path = Int64[]
  topo_sort_recurse( adj, path, collect(1:len) )
end

function topo_sort_recurse( adj::Matrix{Bool}, path::Vector{Int64}, map::Vector{Int64} )
  len = size(adj)[1]
  #println("topo_sort_recurse: len: ",len)
  if len == 0
    return path
  end
  indegree = [ sum(adj[:,i])-1 for i = 1:len ]
  fmin = findmin( indegree )
  if fmin[1] == 0
    r = fmin[2]
    #println("r: ",r,"  map: ",map)
    push!(path,map[r])
    deleteat!( map, r )
    #println("path: ",path,"  map: ",map)
    new_adj = adj[ 1:end .!= r, 1:end .!= r ]  # Remove row r and column r from adj
    topo_sort_recurse( new_adj, path, map )
  end  
  path
end    

# Not used
function remove_vertex( adj::Matrix{Bool}, r::Int64 )
  new_adj = adj[ 1:end .!= r, 1:end .!= r ]
end

# Find all topological sort paths for the adjacency matrix adj
# Returns a double:  (path_list, map_list) where map_list should be a collection of empty lists
function topo_sort_all( adj::Matrix{Bool} )
  len = size(adj)[1]
  topo_sort_recurse_all( adj, Int64[], collect(1:len) )
end

function topo_sort_recurse_all( adj::Matrix{Bool}, path::Vector{Int64}, map::Vector{Int64} )
  len = size(adj)[1]
  #println("topo_sort_recurse_all: len: ",len,"  path: ",path,"  map: ",map)
  if len == 0
    return ([path],[map])
  end
  indegree = [ sum(adj[:,i])-1 for i = 1:len ]
  path_list = Vector{Int64}[]
  map_list = Vector{Int64}[]
  fmin_list = findall( x->x==0, indegree )
  if indegree[fmin_list[1]] == 0
    for r in fmin_list
      #r = fmin_list[1]
      cur_path = deepcopy(path)
      cur_map = deepcopy(map)
      #println("r: ",r,"  map: ",map)
      push!(cur_path,map[r])
      deleteat!( cur_map, r )
      #println("cur_path: ",cur_path,"  cur_map: ",cur_map)
      new_adj = adj[ 1:end .!= r, 1:end .!= r ]
      (plist,mlist) = topo_sort_recurse_all( new_adj, cur_path, cur_map)
      append!(path_list,plist)
      append!(map_list,mlist)
    end
  end
  (path_list,map_list)
end

function topo_sort_check_list( adj::Matrix{Bool}, path_list::Vector{Vector{Int64}}=topo_sort_all(adj)[1] )
  result_list = map( path->topo_sort_check(adj,path), path_list )
  println("result_list: ",result_list)
  if length(findall(x->typeof(x)<:Vector, result_list )) > 0
    return "fail"
  else
    return "pass"
  end
end
  
function topo_sort_check( adj::Matrix{Bool}, path::Vector{Int64}=topo_sort(adj) )
  len = size(adj)[1]
  #path = topo_sort( adj )
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
      
function topo_sorts_phenos( p::Parameters, funcs::Vector{Func}, phlist::GoalList, nreps::Int64; csvfile::String="" )
  max_steps = 50_000
  max_tries = 1
  df = DataFrame( :pheno=> Vector{MyInt}[], :num_topo_sorts=>Float64[] )
  #ngenos = Int128( count_circuits_ch( p, funcs ) )
  for ph in phlist
    nsorts = 0
    count = 0
    for i = 1:nreps
      (ch,step) = neutral_evolution( random_chromosome( p, funcs ), funcs, ph, max_steps)
      ttry = 1
      while ttry < max_tries && step == max_steps
        (ch,step) = neutral_evolution( random_chromosome( p, funcs ), funcs, ph, max_steps)
        ttry += 1
      end
      if step < max_steps
        adj = poset( ch, funcs )
        nsorts += length( topo_sort_all( adj )[1] )
        count += 1
      end
    end
    push!(df, [ ph, (count > 0 ? nsorts/count : 0) ] )
  end

  #=
  if length( phlist ) == 2^2^p.numinputs
    insertcols!(df,3,redundancy_exact( p, funcs ))
  end
  =#
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ",funcs)
      println(f,"# nreps: ",nreps)
      println(f,"# max_steps: ",max_steps)
      println(f,"# max_tries: ",max_tries)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
