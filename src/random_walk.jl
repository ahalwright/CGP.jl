using Statistics, Printf
using DataFrames, CSV, Dates

# Calculates robustness and evolvability for all goal functions by sampling circuits as in Hu et al. (2020).
# The node-adjacency matrix for the network of phenotypes is computed as goal_edge_matrix.
# Diagonal entries are the number of self edges which is the robustness count.
# goal_edge_matrix[g,h] is the number of mutations discovered from goal g to goal h.
# The methodology is to do nprocesses*nwalks random walks each of length steps, and record all of the transitions made.
# If csvfile is specified, writes the dataframe to this file.
function run_random_walks_parallel( nprocesses::Int64, nwalks::Int64, p::Parameters, steps::Int64; 
      csvfile::String="", output_dict::Bool=true )
      #c_list::Vector{Chromosome}=Vector{Chromosome}[] 
  ngoals = 2^(2^p.numinputs)
  if output_dict
    goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
    goal_set_list = pmap( x->run_random_walks( nwalks, p, steps, output_dict=output_dict ), collect(1:nprocesses) )
    println("len: ",length(goal_set_list))
    goal_pair_dict = merge(+,goal_pair_dict,goal_set_list...)
    df = robust_evolvability( goal_pair_dict, p )
  else
    goal_edge_matrix = zeros(Int64,ngoals,ngoals)
    goal_edge_matrix_list = pmap( x->run_random_walks( nwalks, p, steps, output_dict=output_dict ), collect(1:nprocesses) )
    println("len: ",length(goal_edge_matrix_list),"  size: ",size(goal_edge_matrix_list[1]),"  type: ",typeof(goal_edge_matrix_list[1]))
    for gem in goal_edge_matrix_list
      goal_edge_matrix .+= gem
    end
    df = robust_evolvability( goal_edge_matrix )
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# nwalks: ",nwalks)
      println(f,"# steps: ",steps)
      println(f,"# nprocesses: ",nprocesses)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function run_random_walks( nwalks::Int64, p::Parameters, steps::Int64; output_dict::Bool=true) 
  println("run_random_walks with nwalks: ",nwalks,"  steps: ",steps)
  #c_list = Vector{Chromosome}[]
  @assert p.numoutputs == 1
  funcs = default_funcs(p.numinputs)
  ngoals = 2^(2^p.numinputs)
  c_list = [ random_chromosome(p,funcs) for _ = 1:nwalks ]
  if output_dict
    goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
    rw_list = [random_walk( c_list[i], steps, output_dict=output_dict ) for i = 1:nwalks ]
    goal_pair_dict = merge(+,goal_pair_dict,rw_list...)
    return goal_pair_dict
  else
    goal_edge_matrix = zeros(Int64,ngoals,ngoals)
    for i = 1:nwalks
      goal_edge_matrix .+= random_walk( c_list[i], steps )
    end
    return goal_edge_matrix 
  end
end

# Outputs a node-edge adjacency matrix unless output_dict==true, 
#      in which case a dictionary indexed on goal pairs whose values the the count of the goal pair
function random_walk( c::Chromosome, steps::Int64; output_dict::Bool=true )
  funcs = default_funcs(c.params.numinputs)
  ngoals = 2^(2^c.params.numinputs)
  if output_dict
    goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
  else
    goal_edge_matrix = zeros(Int64,ngoals,ngoals)
  end
  goal = output_values(c)[1]
  #println("start goal: ",goal)
  for i = 1:steps 
    (new_c,active) = mutate_chromosome!( c, funcs )
    prev_goal = goal  
    goal = output_values(c)[1]
    #println("i: ",i,"  goal: ",goal)
    if output_dict
      if haskey(goal_pair_dict,(prev_goal,goal))
        goal_pair_dict[(prev_goal,goal)] += 1
      else
        goal_pair_dict[(prev_goal,goal)] = 1
      end
    else
      goal_edge_matrix[Int64(prev_goal)+1,Int64(goal)+1] += 1
    end
  end
  if output_dict
    #println("length(goal_pair_dict: ",length(goal_pair_dict))
    goal_pair_dict
  else
    goal_edge_matrix
  end
end

# Calculates robustness and degree evolvability for each goal and saves these in a dataframe.
function robust_evolvability( goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, p::Parameters )
  ngoals = 2^(2^p.numinputs)
  robustness = zeros(Float64,ngoals)
  evolvability = zeros(Float64,ngoals)
  all_goals = MyInt(0):MyInt(ngoals-1)
  dget(x) = get(goal_pair_dict,x,0)  # get from goal_pair_dict with default value 0
  onez(x) = x > 0 ? 1 : 0   # function that maps anything greater than zero to 1
  for g = MyInt(0):MyInt(ngoals-1)
    robustness[g+1] = dget((g,g))
    evolvability[g+1] = sum( onez(dget((g,h))+dget((h,g))) for h = all_goals ) - onez(robustness[g+1])
  end
  df = DataFrame()
  df.goal = collect(all_goals)
  df.robustness = robustness
  df.evolvability = evolvability
  df
end

# Calculates robustness and degree evolvability for each goal and saves these in a dataframe.
function robust_evolvability( goal_edge_matrix::Array{Int64,2} )
  ngoals = size(goal_edge_matrix)[1]
  triangularize!(goal_edge_matrix)
  robustness = zeros(Float64,ngoals)
  evolvability = zeros(Float64,ngoals)
  for g = 1:ngoals
    robustness[g] = goal_edge_matrix[g,g]
    sum1 = g==ngoals ? 0 : sum( ((goal_edge_matrix[g,h] != 0) ? 1 : 0) for h = (g+1):ngoals )
    pairs1 = g==ngoals ? 0 : [ (g,h) for h = (g+1):ngoals ] 
    sum2 = g==1 ? 0 : sum( ((goal_edge_matrix[h,g] != 0) ? 1 : 0) for h = 1:(g-1))
    pairs2 = g==1 ? 0 : [ ((h,g),goal_edge_matrix[h,g]) for h = 1:(g-1)] 
    #println("g: ",g,"  pairs1: ",pairs1,"  sum1: ",sum1)
    #println("g: ",g,"  pairs2: ",pairs2,"  sum2: ",sum2)
    evolvability[g] = sum1 + sum2
  end
  df = DataFrame()
  df.goal = collect(MyInt(0):MyInt(ngoals-1))
  df.robustness = robustness
  df.evolvability = evolvability
  df
end

# Convert into a upper triangular matrix by adding the below-diagonal entries to the corresponding above-diagonal entry
function triangularize!( M::Array{Int64,2} )
  ngoals = size(M)[1]
  @assert ngoals == size(M)[2]
  for g = 1:(ngoals-1)
    for h = (g+1):ngoals
      M[g,h] += M[h,g]
      M[h,g] = 0
    end
  end   
end

# Convert into a upper triangular matrix by adding the below-diagonal entries to the corresponding above-diagonal entry
function triangularize!( M::Array{Int64,2} )
  ngoals = size(M)[1]
  @assert ngoals == size(M)[2]
  for g = 1:(ngoals-1)
    for h = (g+1):ngoals
      M[g,h] += M[h,g]
      M[h,g] = 0
    end
  end
end

function matrix_to_dict( M::Array{Int64,2} )
  goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
  for g = MyInt(0):(size(M)[1]-1)
    for h = MyInt(0):(size(M)[2]-1)
      if M[g+1,h+1] != 0
        goal_pair_dict[(g,h)] = M[g+1,h+1]
      end
    end
  end
  goal_pair_dict
end 

function dict_to_matrix( d::Dict{Tuple{MyInt,MyInt}}, size::Int64 )
  M = zeros(Int64,size,size)
  for k in keys(d)
    (g,h) = k
    M[g+1,h+1] = d[k]
  end
  M
end
