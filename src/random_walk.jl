# Used in GECCO paper
using Statistics, Printf
using DataFrames, CSV, Dates
# This file does not include the "export" declarations needed to be included in CGP.jl.
# Thus, it should be directly included.  Sample runs in ../data/9_3_21/.

# Calculates robustness and evolvability for all phenotypes by sampling circuits as in Hu et al. (2020).
# The node-adjacency matrix for the network of phenotypes is computed as goal_edge_matrix.
# Diagonal entries are the number of self edges which is the robustness count.
# goal_edge_matrix[g,h] is the number of mutations discovered from goal g to goal h.
# The methodology is to do nprocesses*nwalks random walks each of length steps, and record all of the transitions made.
# gl is a list of MyInts (for single-output goals) that robust_evolvability() uses to create the output dataframe which 
#   includes robustness, d_evolvability, and s_evolvability.
# If csvfile is specified, writes the dataframe to this file.
function run_random_walks_parallel( nprocesses::Int64, nwalks::Int64, gl::Vector{MyInt}, p::Parameters, steps::Int64; 
      csvfile::String="", output_dict::Bool=true, save_complex::Bool=false, use_lincircuit::Bool=false )
  if save_complex && !output_dict
    error("The save_complex=true and output_dict=false options are not compatible in run_random_walks_parallel().")
  end
  ngoals = 2^(2^p.numinputs)
  addvalues(x::Tuple{Int64,Float64},y::Tuple{Int64,Float64}) = (x[1]+y[1],(x[2]+y[2])/2.0)
  if output_dict
    goal_pair_dict = save_complex ? Dict{Tuple{MyInt,MyInt},Tuple{Int64,Float64}}() : Dict{Tuple{MyInt,MyInt},Int64}()
    #goal_pair_dict = Dict{Tuple{MyInt,MyInt},Int64}()
    dict_list = pmap( x->run_random_walks( nwalks, p, steps, output_dict=output_dict, save_complex=save_complex, use_lincircuit=use_lincircuit ), 
        collect(1:nprocesses) )
    #println("length(dict_list): ",length(dict_list))
    #println("dict_list[1]: ",dict_list[1])
    #println("typeof(goal_pair_dict): ",typeof(goal_pair_dict))
    #for gp in dict_list
      goal_pair_dict = save_complex ? merge(addvalues,goal_pair_dict,dict_list...) : merge(+,goal_pair_dict,dict_list...)
      #goal_pair_dict = merge(+, goal_pair_dict, gp )
    # end
    #return goal_pair_dict
    df = robust_evolvability( goal_pair_dict, gl, p, save_complex )  # The version of robust_evolvability() called depends on the type of goal_pair_dict
  else
    goal_edge_matrix = zeros(Int64,ngoals,ngoals)
    goal_edge_matrix_list = pmap( x->run_random_walks( nwalks, p, steps, output_dict=output_dict, use_lincircuit=use_lincircuit ), 
        collect(1:nprocesses) )
    println("len: ",length(goal_edge_matrix_list),"  size: ",size(goal_edge_matrix_list[1]),"  type: ",typeof(goal_edge_matrix_list[1]))
    for gem in goal_edge_matrix_list
      goal_edge_matrix .+= gem
    end
    #println(goal_edge_matrix)
    df = robust_evolvability( goal_edge_matrix, gl )
  end
  if length(csvfile) > 0
    println("csvfile: ",csvfile)
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# nwalks: ",nwalks)
      println(f,"# steps: ",steps)
      println(f,"# nprocesses: ",nprocesses)
      println(f,"# output_dict: ",output_dict)
      println(f,"# use_lincircuit: ",use_lincircuit)
      println(f,"# save_complex: ",save_complex)
      CSV.write( f, df, append=true, writeheader=true )
      println("csvfile written")
    end
  end
  df
  #goal_edge_matrix
end

function run_random_walks( nwalks::Int64, p::Parameters, steps::Int64; use_lincircuit::Bool=false, output_dict::Bool=false, save_complex::Bool=false )
  @assert p.numoutputs == 1
  funcs = default_funcs(p.numinputs)
  addvalues(x::Tuple{Int64,Float64},y::Tuple{Int64,Float64}) = (x[1]+y[1],(x[2]+y[2])/2.0)
  ngoals = 2^(2^p.numinputs)
  if use_lincircuit
    c_list = [ rand_lcircuit(p,funcs) for _ = 1:nwalks ]
  else
    c_list = [ random_chromosome(p,funcs) for _ = 1:nwalks ]
  end
  if output_dict
    goal_pair_dict = save_complex ? Dict{Tuple{MyInt,MyInt},Tuple{Int64,Float64}}() : Dict{Tuple{MyInt,MyInt},Int64}()
    for i = 1:nwalks
      random_walk!( goal_pair_dict, c_list[i], steps )
    end
    return goal_pair_dict
  else
    goal_edge_matrix = zeros(Int64,ngoals,ngoals)
    for i = 1:nwalks
      random_walk!( goal_edge_matrix, c_list[i], steps )
      #println(sum(goal_edge_matrix))
    end
    return goal_edge_matrix 
  end
end

# Does one random walk starting at c
# Outputs a node-edge adjacency matrix unless output_dict==true, 
#      in which case a dictionary indexed on goal pairs whose values the the count of the goal pair
function random_walk!( dict_mat::Union{Dict{Tuple{MyInt,MyInt},Tuple{Int64,Float64}}, Dict{Tuple{MyInt,MyInt},Int64},Matrix{Int64}},
    c::Union{Chromosome,LinCircuit}, steps::Int64 )
  funcs = default_funcs(c.params.numinputs)
  ngoals = 2^(2^c.params.numinputs)
  #addvalues(x,y) = (x[1]+y[1],(x[2]+y[2])/2.0)
  addvalues(x::Tuple{Int64,Float64},y::Tuple{Int64,Float64}) = (x[1]+y[1],(x[2]+y[2])/2.0)
  output_dict = false  # establish scope
  save_complex = false  # establish scope
  if typeof(dict_mat) == Matrix{Int64}
    goal_edge_matrix = dict_mat
    output_dict = false
  elseif typeof(dict_mat) == Dict{Tuple{MyInt,MyInt},Tuple{Int64,Float64}}
    goal_pair_dict = dict_mat
    output_dict = true
    save_complex = true
  elseif typeof(dict_mat) == Dict{Tuple{MyInt,MyInt},Int64}
    goal_pair_dict = dict_mat
    output_dict = true
    save_complex = false
  end
  #println("output_dict: ",output_dict,"  save_complex: ",save_complex)
  cmplx = 0.0
  goal = output_values(c)[1]
  #println("start goal: ",goal)
  for i = 1:steps 
    prev_goal = goal
    #println("i: ",i,"  prev_goal: ",prev_goal,"  goal: ",goal)
    if typeof(c) == Chromosome
      cmplx = save_complex ? complexity5(c) : 0.0 
      mutate_chromosome!( c, funcs )
      #print_circuit(c)
    elseif typeof(c) == LinCircuit
      cmplx = save_complex ? lincomplexity( c, funcs )  : 0
      mutate_circuit!( c, funcs )
      #print_circuit(c, funcs )
    end
    goal = output_values(c)[1]
    if output_dict 
      value = save_complex ? (1,cmplx) : 1
      if haskey(goal_pair_dict,(prev_goal,goal))
        prev_value = goal_pair_dict[(prev_goal,goal)]
        goal_pair_dict[(prev_goal,goal)] = save_complex ?  addvalues(prev_value,value) : prev_value+value
      else
        goal_pair_dict[(prev_goal,goal)] = value
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
function robust_evolvability( goal_pair_dict::Dict{Tuple{MyInt,MyInt},Tuple{Int64,Float64}}, gl::Vector{MyInt}, p::Parameters, save_complex::Bool )
  ngoals = 2^(2^p.numinputs)
  allgoals = MyInt(0):MyInt(ngoals-1)
  println("allgoals: ",allgoals,"  ngoals: ",ngoals)
  #re_list = pmap( g->robust_evo(g, goal_pair_dict, ngoals ), gl )
  re_list = map( g->robust_evo(g, goal_pair_dict, ngoals ), gl )
  df = DataFrame()
  df.goal = gl
  df.frequency= [ re[1] for re in re_list ]
  df.robustness = [ re[2] for re in re_list ]
  df.s_evolvability = [ re[3] for re in re_list ]  
  df.d_evolvability = [ re[4] for re in re_list ]  
  if save_complex
    df.complexity= [ re[5] for re in re_list ]  
  end
  df
end

# Helper function used in robust_evolvability()
# Compute output dataframe
function robust_evo( g::MyInt, goal_pair_dict::Dict{Tuple{MyInt,MyInt},Tuple{Int64,Float64}}, ngoals::Int64 )
  allgoals = MyInt(0):MyInt(ngoals-1)
  dget(x) = get(goal_pair_dict,x,(0,0.0))
  dget10(x) = get(goal_pair_dict,x,(0,0.0))[1]>0 ? 1 : 0
  #println("dget((g,g)): ",dget((g,g)),"  dget((g,g))[2]: ",dget((g,g))[2])
  sum_gh = sum( dget((g,h))[1] for h in allgoals )
  #sum_hg = sum( dget((h,g))[1] for h in allgoals )
  #for h in allgoals println("  dget((g,h)): ",dget((g,h))) end
  sum_gh_c = sum( dget((g,h))[1]*dget((g,h))[2] for h in allgoals )
  gg = dget((g,g))[1]
  frequency = sum_gh
  robustness = gg/frequency    
  s_evolvability = sum_gh - gg  # count of mutations that take g to a different phenotype
  d_evolvability = sum(dget10((g,h)) for h in allgoals) - ((gg>0) ? 1 : 0)  #count of number of phenotypes by mutation of g
  complexity = sum_gh_c/frequency
  (frequency,robustness,s_evolvability,d_evolvability,complexity)
end

# Note that gl is a list of MyInts, not a list of goals 
function robust_evolvability( goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, gl::Vector{MyInt}, p::Parameters, save_complex::Bool=false )
  ngoals = 2^(2^p.numinputs)
  allgoals = MyInt(0):MyInt(ngoals-1)
  println("allgoals: ",allgoals)
  #re_list = pmap( g->robust_evo(g, goal_pair_dict, ngoals ), gl )
  re_list = map( g->robust_evo(g, goal_pair_dict, ngoals ), gl )
  df = DataFrame()
  df.goal = gl
  df.frequency= [ re[1] for re in re_list ]
  df.robustness = [ re[2] for re in re_list ]
  df.s_evolvability = [ re[3] for re in re_list ]  
  df.d_evolvability = [ re[4] for re in re_list ]  
  df
end

# Helper function for robust_evolvability()
# Corrected 1/2/21 and 1/3/21
function robust_evo( g::MyInt, goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, ngoals::Int64 )
  allgoals = MyInt(0):MyInt(ngoals-1)
  dget(x) = get(goal_pair_dict,x,0)
  dget10(x) = get(goal_pair_dict,x,0)>0 ? 1 : 0
  sum_gh = sum( dget((g,h)) for h in allgoals )
  #sum_hg = sum( dget((h,g)) for h in allgoals )
  gg = dget((g,g))
  frequency = sum_gh
  robustness = gg/frequency    
  s_evolvability = sum_gh - gg  # count of mutations that take g to a different phenotype
  d_evolvability = sum(dget10((g,h)) for h in allgoals) - ((gg>0) ? 1 : 0)  #count of number of phenotypes by mutation of g
  (frequency,robustness,s_evolvability,d_evolvability)
end

# Calculates robustness and degree evolvability for each goal and saves these in a dataframe.
# Robustness 
function robust_evolvability( goal_edge_matrix::Array{Int64,2}, gl::Vector{MyInt} )
  ngoals = size(goal_edge_matrix)[1]
  #triangularize!(goal_edge_matrix)
  frequency = zeros(Int64,ngoals)
  robustness = zeros(Float64,ngoals)
  s_evolvability = zeros(Int64,ngoals)
  d_evolvability = zeros(Int64,ngoals)
  for g in map(x->x+1,gl)   # convert to 1-based indexing
    frequency[g] = sum( goal_edge_matrix[g,h]  for h = 1:ngoals )
    robustness[g] = goal_edge_matrix[g,g]/frequency[g]
    s_evolvability[g] = frequency[g] - goal_edge_matrix[g,g]
    d_evolvability[g] = sum( ((goal_edge_matrix[g,h] != 0) ? 1 : 0) for h = 1:ngoals ) - ((goal_edge_matrix[g,g]>0) ? 1 : 0)
  end
  df = DataFrame()
  df.goal = collect(MyInt(0):MyInt(ngoals-1))
  df.frequency = frequency
  df.robustness = robustness
  df.s_evolvability = s_evolvability
  df.d_evolvability = d_evolvability
  df
end

# Convert into a upper triangular matrix by adding the below-diagonal entries to the corresponding above-diagonal entry
# Never used.
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
  for g = 0:(size(M)[1]-1)
    for h = 0:(size(M)[2]-1)
      if M[g+1,h+1] != 0
        goal_pair_dict[(g,h)] = M[g+1,h+1]
      end
    end
  end
  goal_pair_dict
end

function dict_to_matrix( goal_pair_dict::Dict{Tuple{MyInt,MyInt},Int64}, p::Parameters )
  M = zeros(Int64,2^2^p.numinputs,2^2^p.numinputs)
  for (g,h) in keys(goal_pair_dict)
    M[g+1,h+1] = goal_pair_dict[(g,h)] 
  end
  M
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
  robustness = zeros(Float64,length(gl))
  s_evolvability = zeros(Int64,length(gl))
  d_evolvability = zeros(Int64,length(gl))
  allgoals = MyInt(0):MyInt(2^2^p.numinputs-1)
  i = 1
  for g in gl
    #frequency[i] = sum( get( goal_pair_dict, (g,h), 0 ) + get( goal_pair_dict, (h,g), 0 ) for h in allgoals ) 
    frequency[i] = sum( dget((g,h)) for h in allgoals ) 
    #robustness[i] = get( goal_pair_dict, (g,g), 0 )
    robustness[i] = dget((g,g))/frequency[i]
    s_evolvability[i] = sum( dget((g,h)) for h in allgoals ) - dget((g,g))
    d_evolvability[i] = sum( dget10((g,h)) for h in allgoals ) - (dget((g,g))>0 ? 1 : 0)
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
