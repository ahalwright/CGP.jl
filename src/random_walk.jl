using Statistics, Printf
using DataFrames, CSV, Dates

# Calculates robustness and evolvability for all goal functions by sampling circuits as in Hu et al. (2020).
# The node-adjacency matrix for the network of phenotypes is computed as goal_edge_matrix.
# Diagonal entries are the number of self edges which is the robustness count.
# goal_edge_matrix[g,h] is the number of mutations discovered from goal g to goal h.
# The methodology is to do nprocesses*nwalks random walks each of length steps, and record all of the transitions made.
# If csvfile is specified, writes the dataframe to this file.
function run_random_walks_parallel( nprocesses::Int64, nwalks::Int64, gl::Vector{MyInt}, p::Parameters, steps::Int64; 
      csvfile::String="", output_dict::Bool=true )
      #c_list::Vector{Chromosome}=Vector{Chromosome}[] 
  ngoals = 2^(2^p.numinputs)
  if output_dict
    goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
    goal_set_list = pmap( x->run_random_walks( nwalks, p, steps, output_dict=output_dict ), collect(1:nprocesses) )
    println("length(goal_set_list): ",length(goal_set_list))
    for gp in goal_set_list
      goal_pair_dict = merge(+, goal_pair_dict, gp )
    end
    #return goal_pair_dict
    df = robust_evolvability( goal_pair_dict, gl, p )
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

function run_random_walks( nwalks::Int64, p::Parameters, steps::Int64; output_dict::Bool=true )
      #output_dict::Bool=true, c_list::Vector{Chromosome}=Vector{Chromosome}[] )
  #println("run_random_walks with nwalks: ",nwalks,"  steps: ",steps)
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
function random_walk( c::Chromosome, steps::Int64; output_dict::Bool=true, mutate_locs::Vector{Int64}=Int64[] )
  #mutate_locs = Int64[]
  funcs = default_funcs(c.params.numinputs)
  ngoals = 2^(2^c.params.numinputs)
  loc_index = 1
  if length(mutate_locs) == 0
    #nlocs = length(mutate_all(deepcopy(c),funcs))
    nlocs = num_mutate_locations( c, funcs )
    mutate_locs = rand(1:nlocs,steps)
  end
  if output_dict
    goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
  else
    goal_edge_matrix = zeros(Int64,ngoals,ngoals)
  end
  goal = output_values(c)[1]
  #println("start goal: ",goal)
  for i = 1:steps 
    #println("loc_index: ",loc_index,"  mutate_locs[loc_index]: ",mutate_locs[loc_index])
    mutate_chromosome!( c, funcs, mutate_locs[loc_index])
    #print_circuit(c)
    #print_circuit(new_c)
    loc_index += 1
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
    #println("length(goal_pair_dict): ",length(goal_pair_dict))
    goal_pair_dict
  else
    goal_edge_matrix
  end
end

function robust_evolvability( goal_pair_list::Array{Pair{Tuple{UInt16,UInt16},Int64},1}, p::Parameters)
  robust_evolvability( pairs_to_dict( goal_pair_list ), p )
end

# 
# Note that gl is a list of MyInts, not a list of goals 
function robust_evolvability( goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, gl::Vector{MyInt}, p::Parameters )
  ngoals = 2^(2^p.numinputs)
  allgoals = MyInt(0):MyInt(ngoals-1)
  re_list = pmap( g->robust_evo(g, goal_pair_dict, ngoals ), gl )
  #re_list = map( g->robust_evo(g, goal_pair_dict, ngoals ), gl )
  df = DataFrame()
  df.goal = gl
  df.frequency= [ re[1] for re in re_list ]
  df.robustness = [ re[2] for re in re_list ]
  df.s_evolvability = [ re[3] for re in re_list ]  
  df.d_evolvability = [ re[4] for re in re_list ]  
  df
end

# Helper function used in pmap() in robust_evolvability()
function robust_evo( g::MyInt, goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, ngoals::Int64 )
  allgoals = MyInt(0):MyInt(ngoals-1)
  dget(x) = get(goal_pair_dict,x,0)
  dget10(x) = get(goal_pair_dict,x,0)>0 ? 1 : 0
  frequency = sum( dget((g,h)) + dget((h,g)) for h in allgoals )
  robustness = get(goal_pair_dict,(g,g),0)
  s_evolvability = sum( dget((g,h)) + dget((h,g)) for h = allgoals ) - 2*robustness
  d_evolvability = sum( dget10((g,h)) + dget10((h,g)) for h = allgoals ) - 2*(robustness>0 ? 1 : 0)
  (frequency,robustness,s_evolvability,d_evolvability)
end

# Calculates robustness and degree evolvability for each goal and saves these in a dataframe.
# Robustness 
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

function matrix_to_dict( M::Array{Int64,2} )
  goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
  for g = 1:size(M)[1]
    for h = 1:size(M)[2]
      if M[g,h] != 0
        goal_pair_dict[(g,h)] = M[g,h]
      end
    end
  end
  goal_pair_dict
end

function pairs_to_dict( pairs::Array{Pair{Tuple{MyInt,MyInt},Int64},1})
  println("pairs to dict")
  pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
  for p in pairs
    pair_dict[p[1]] = p[2]
  end
  pair_dict
end 

# Run robust_evolvability and then write to csvfile if it is given
# Ran for over 24 hours on a large goal_pair_list.
function robust_evolvability_to_df( goal_pair_list::Array{Pair{Tuple{UInt16,UInt16},Int64},1}, p::Parameters; 
      csvfile::String="" ) 
  df = robust_evolvability( pairs_to_dict( goal_pair_list ), p )
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      #println(f,"# nwalks: ",nwalks)
      #println(f,"# steps: ",steps)
      #println(f,"# nprocesses: ",nprocesses)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Computes frequency, robustess, strength_evolvability, degree_evolvabilty of a list of MyInts (interpreted as goals)

# Note that gl is a list of MyInts, not a list of goals 
function dict_to_evolvability( goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, gl::Vector{MyInt}, p::Parameters )
  @assert p.numoutputs == 1
  dget(x) = get(goal_pair_dict,x,0)
  dget10(x) = get(goal_pair_dict,x,0)>0 ? 1 : 0
  frequency = zeros(Int64,length(gl))
  robustness = zeros(Int64,length(gl))
  s_evolvability = zeros(Int64,length(gl))
  d_evolvability = zeros(Int64,length(gl))
  allgoals = MyInt(0):MyInt(2^2^p.numinputs-1)
  i = 1
  for g in gl
    frequency[i] = sum( get( goal_pair_dict, (g,h), 0 ) + get( goal_pair_dict, (h,g), 0 ) for h in allgoals ) 
    robustness[i] = get( goal_pair_dict, (g,g), 0 )
    s_evolvability[i] = sum( dget((g,h)) + dget((h,g)) for h in allgoals ) - 2*robustness[i]
    d_evolvability[i] = sum( dget10((g,h)) + dget10((h,g)) for h in allgoals ) - 2*(robustness[i]>0 ? 1 : 0)
    i += 1
  end
  df = DataFrame()
  df.goal = gl
  df.frequency = frequency
  df.robustness = robustness
  df.s_evolvability = s_evolvability
  df.d_evolvability = d_evolvability
  df
end
