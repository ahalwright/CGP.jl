# Includes explore_complexity()
using DataFrames, Statistics, Dates, CSV, Distributed, Printf, Random, Plots
export evolve_complexity, run_complexity_list, complexity_list
export run_explore_complexity, explore_complexity
export average_dataframes, extend_list_by_dups, goal_complexity_frequency_dataframe
export filter_goallist_by_complexity, complexity_freq_scatter_plot
export bin_value, bin_data, bin_counts
export kolmogorov_complexity, run_kolmogorov_complexity, redund_vs_k_complexity_plot
export run_k_complexity_mutate_all, kcomp_summary_dataframe, run_redundancy_mutate_all 
export redundancy_dict, kolmogorov_complexity_dict, testKall
##gti, gti!, 
# Distribution of complexities for circuits evolving into a given goal
# Results in data/10_31
function evolve_complexity( g::Goal, p::Parameters, num_circuits::Int64, maxsteps::Int64 )
  n_repeats = 20
  funcs = default_funcs(p.numinputs)
  total_steps = 0
  complexity_list = Float64[]
    c = random_chromosome(p,funcs)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [g], funcs, maxsteps )  
    total_steps += step
    count = 0
    while step == maxsteps && count < n_repeats
      (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [g], funcs, maxsteps )
      total_steps += step
      #println("total_steps: ",total_steps)
      count += 1
    end   
    if step == maxsteps || count == n_repeats
      println("step == maxsteps in function evolve_to_random_goals()")  
    end   
    cmplx = complexity5(c)
  #println("(cmplx,total_steps): ",(cmplx,total_steps))
  (cmplx,total_steps)
end

function run_complexity_list( ngoals::Int64, p::Parameters, num_circuits::Int64, maxsteps::Int64; csvfile::String="" )
  goallist = randgoallist( ngoals, p.numinputs, p.numoutputs )
  run_complexity_list(goallist, p, num_circuits,maxsteps, csvfile=csvfile )
end

function run_complexity_list( goallist::GoalList, p::Parameters, num_circuits::Int64, maxsteps::Int64; csvfile::String="" )
  df = DataFrame()
  df.goal = Goal[]
  df.steps = Float64[]
  df.mean = Float64[]
  df.sddev = Float64[]
  df.max = Float64[]
  df.min = Float64[]
  df.q10 = Float64[]
  df.q90 = Float64[]
  df.q95 = Float64[]
  df.q90 = Float64[]
  df.q95 = Float64[]
  df.q99 = Float64[]
  extended_goallist = vcat([ fill(g,num_circuits) for g in goallist ]... )
  #println("length(extended_goallist): ",length(extended_goallist),"  extended_goallist: ",extended_goallist)
  complexity_lists = pmap( g->evolve_complexity( g, p, num_circuits, maxsteps ), extended_goallist )
  #println("complexity_lists: ",complexity_lists)
  ngoals = length(goallist)
  for i = 1:ngoals
    complexity_list = Float64[]
    steps_list = Int64[]
    for j = 1:num_circuits
      #println("(i,j): ",(i,j),"  ind: ",(i-1)*num_circuits+j)
      push!(complexity_list,complexity_lists[(i-1)*num_circuits+j][1])
      push!(steps_list,complexity_lists[(i-1)*num_circuits+j][2])
    end
    #println("steps_list: ",steps_list)
    row_tuple = 
      (
        goallist[i],
        mean(steps_list),
        mean(complexity_list),
        std(complexity_list),
        maximum(complexity_list),
        minimum(complexity_list),
        quantile(complexity_list,0.10),
        quantile(complexity_list,0.90),
        quantile(complexity_list,0.95),
        quantile(complexity_list,0.99)
      )
    #println("row_tuple: ",row_tuple)
    push!(df,row_tuple)
  end
  hostname = readchomp(`hostname`)
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# maxsteps: ",maxsteps)
      println(f,"# num_circuits: ",num_circuits)
      CSV.write( f, df, append=true, writeheader=true )
    end
  else
    println("# date and time: ",Dates.now())
    println("# host: ",hostname," with ",nprocs()-1,"  processes: " )
    println("# funcs: ", Main.CGP.default_funcs(p.numinputs))
    print_parameters(p,comment=true)
    println("# maxsteps: ",maxsteps)
      println("# num_circuits: ",num_circuits)
  end
  df
end

# Run the next version of run_explore_complexity() nreps time parallelized by pmap().
# This version generates the goallist.
function run_explore_complexity( nreps::Int64, nruns::Int64, p::Parameters, ngoals::Int64, num_circuits::Int64, max_ev_steps::Int64;
      csvfile::String="", insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  goallist = randgoallist( ngoals, p.numinputs,p.numoutputs)
  run_explore_complexity( nreps, nruns, p, goallist, num_circuits, max_ev_steps, csvfile=csvfile )
end

# Run the next version of run_explore_complexity() nreps time parallelized by pmap().
# This version takes goallist as an argument.
function run_explore_complexity( nreps::Int64, nruns::Int64, p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64;
      csvfile::String="", insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  ngoals = length(goallist)
  result_dfs = DataFrame[]
  #result_dfs = pmap(x->run_explore_complexity( nruns, p, goallist, num_circuits, max_ev_steps,
  #      insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob), collect(1:nreps))
  result_dfs = map(x->run_explore_complexity( nruns, p, goallist, num_circuits, max_ev_steps, 
        insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob), collect(1:nreps))
  println("length(result_dfs): ",length(result_dfs))
  df = average_dataframes( result_dfs )
  insertcols!(df, 1, :generation=>collect(1:size(df)[1]))
  cumm_unique_goals = zeros(Float64,size(df)[1])
  cumm_unique_goals[1] = df.unique_goals[1]
  for i = 2:size(df)[1]
    cumm_unique_goals[i] = df.unique_goals[i] + cumm_unique_goals[i-1]
  end
  insertcols!(df, size(df)[2]+1, :cumm_unique_goals=>cumm_unique_goals )
  hostname = readchomp(`hostname`)
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# nreps: ",nreps)
      println(f,"# nruns: ",nruns)
      println(f,"# num_circuits: ",num_circuits)
      println(f,"# max_ev_steps: ",max_ev_steps)
      println(f,"# ngoals: ",ngoals)
      if insert_gate_prob > 0.0
        println(f,"# insert_gate_prob: ",insert_gate_prob)
      end 
      if delete_gate_prob > 0.0
        println(f,"# delete_gate_prob: ",delete_gate_prob)
      end 
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  println("# date and time: ",Dates.now())
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )
  println("# funcs: ", Main.CGP.default_funcs(p.numinputs))
  print_parameters(p,comment=true)
  println("# nreps: ",nreps)
  println("# nruns: ",nruns)
  println("# num_circuits: ",num_circuits)
  println("# max_ev_steps: ",max_ev_steps)
  println("# ngoals: ",ngoals)
  if insert_gate_prob > 0.0
    println("# insert_gate_prob: ",insert_gate_prob)
  end 
  if delete_gate_prob > 0.0
    println("# delete_gate_prob: ",delete_gate_prob)
  end 
  df
end

# Run explore_complexity() nruns times with the discovered circuit_list from one run becoming the starting circuit_list for the next.
# Also, discoverd goals are removed from goallist on subsequent runs.
function run_explore_complexity( nruns::Int64, p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64;
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  max_complexity_tries = 40
  df = DataFrame()
  done = false
  explore_complexity_tries = 0
  while !done && explore_complexity_tries <= max_complexity_tries
    explore_complexity_tries += 1
    #println("while !done loop.")
    df = DataFrame()
    df.num_gates = Float64[]
    df.levelsback = Float64[]
    #=  Too complicated for now.  5/19/22
    if insert_gate_prob > 0.0
      df.insert_gate_prob = Float64[]
    end 
    if delete_gate_prob > 0.0
      df.delte_gate_prob = Float64[]
    end 
    =#
    df.steps = Int64[]
    df.unique_goals = Int64[]
    df.mean_complexity = Float64[]
    df.max_complexity = Float64[]
    (circuit_list,goals_found,row_tuple) = explore_complexity( p, goallist, num_circuits, max_ev_steps,
        insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
    push!(df,row_tuple)
    reduced_goallist = goallist
    for i = 2:(nruns-1)
      extended_circuit_list = extend_list_by_dups( circuit_list, num_circuits )
      reduced_goallist = setdiff(reduced_goallist,goals_found)
      #println("run i: ",i,"  length(circuit_list): ",length(circuit_list),"  length(reduced_goallist): ",length(reduced_goallist))
      (circuit_list,goals_found,row_tuple) = explore_complexity( p, reduced_goallist, extended_circuit_list, max_ev_steps,
        insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
      if row_tuple[1] == 0
        println("i: ",i," break")
        break # Try again
      end
      push!(df,row_tuple)
    end
    extended_circuit_list = extend_list_by_dups( circuit_list, num_circuits )
    reduced_goallist = setdiff(reduced_goallist,goals_found)
    (circuit_list,goals_found,row_tuple) = explore_complexity( p, reduced_goallist, extended_circuit_list, max_ev_steps,
        insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
    #println("length(redued_goallist): ",length(reduced_goallist),"  row_tuple: ",row_tuple)
    push!(df,row_tuple)
    if row_tuple[3] > 0
      done = true
    end
  end
  println("explore complexity_tries: ",explore_complexity_tries)
  df
end

# This version generates circuit list by calling random_chromosome() to generate each circuit.
function explore_complexity( p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64;
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  funcs = default_funcs(p.numinputs)
  circuit_list = [ random_chromosome(p,funcs) for _=1:num_circuits ]
  explore_complexity( p, goallist, circuit_list, max_ev_steps, 
      insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
end
  
# This is the explore_complexity() function that runs one generation starting with circuits in circuit_list.
# max_ev_steps  is the maximum number of steps in the evolution of a circuit to any of the goals in goalllist.
# Returns a triple (circuits_list, goals_found, row_tuple) where:
#   circuits_list is the list of evolved circuits that output one of the goals
#   goals_found is the list of goals that are found by evolution
#   row_tuple is the row of statistics that can be pushed onto a dataframe
function explore_complexity( p::Parameters, goallist::GoalList, circuit_list::Vector{Chromosome}, max_ev_steps::Int64;
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  #println("explore_complexity  insert_gate_prob: ",insert_gate_prob)
  funcs = default_funcs(p.numinputs)
  goals_found = Goal[]
  steps_list = Int64[]
  src_cmplx_list = Float64[]
  dest_cmplx_list = Float64[]
  circuits_list = Chromosome[]
  num_gates_sum = 0
  levelsback_sum = 0
  for c in circuit_list
    (new_c,step,worse,same,better,output,matched_goals,matched_goals_list) = 
        mut_evolve( c, goallist, funcs, max_ev_steps, print_steps=false, insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
    if step < max_ev_steps
      #print("i: ",i,"  steps: ",step)
      #@printf("  goal: 0x%04x\n",matched_goals_list[1][1][1])
      num_gates_sum += new_c.params.numinteriors
      levelsback_sum += new_c.params.numlevelsback
      push!(goals_found, [matched_goals_list[1][1][1]])
      push!(steps_list,step)
      push!(src_cmplx_list, complexity5(c))
      push!(dest_cmplx_list, complexity5(new_c))
      push!(circuits_list, new_c )
    end
  end
  if length(steps_list) > 0
    num_gates_avg = num_gates_sum/length(steps_list)
    levelsback_avg = levelsback_sum/length(steps_list)
    row_tuple = (num_gates_avg, levelsback_avg, length(steps_list), length(unique(goals_found)),mean(dest_cmplx_list), maximum(dest_cmplx_list) )
  else
    row_tuple = (p.numinteriors, p.numlevelsback,0,0,0.0,0.0)
  end
  #println("length(circuits_list): ",length(circuits_list))
  return (circuits_list, goals_found, row_tuple)
end

# Example:  pops_to_tbl([2,2,1,3],[1,1,1,3]) =
#  2×3 Array{Float64,2}:
#   0.25  0.125  0.125
#   0.0   0.375  0.125
function pops_to_tbl( P::Vector{CGP.FPopulation} )
  lengths = map(length,P)
  @assert( all(lengths .==  lengths[1]))   # Check that esch population has the same length
  sums = map(sum,P)
  try
    @assert( all(sums .≈ 1.0 ))   # Check that each population has approximate sum 1.0
  catch
    println("assertion error pops_to_tbl:  sums: ",sums)
  end
  result = zeros(Float64,length(P),length(P[1]))
  for i = 1:length(P)
    result[i,:] = P[i]/length(P)
  end
  result
end  

# Averages a list of dataframes which should all be the same size with the same names.
# All columns should be numerical.  All columns of the result are Float64.
function average_dataframes( dataframes_list::Vector{DataFrame} )
  if length(dataframes_list) == 0
    return nothing
  end
  df = DataFrame()
  i = 1
  for nm in names(dataframes_list[1])
    insertcols!(df,i,nm=>zeros(Float64,size(dataframes_list[1])[1]))
    i += 1
  end
  #println("df: ",df)
  for d in dataframes_list
    @assert size(d) == size(dataframes_list[1])
    @assert names(d) == names(dataframes_list[1])
    #println("  d: ",d)
    for i = 1:size(dataframes_list[1])[2] 
      df[:,i] = df[:,i] + d[:,i]
    end
  end
  for i = 1:size(dataframes_list[1])[2] 
    df[:,i] = df[:,i]/length(dataframes_list)
  end
  df
end

# Extend vector v to new_length by concatenating v with itself.
function extend_list_by_dups( v::AbstractVector, new_length::Int64 )
  len = length(v)
  if len == 0
    return v
  end
  sum_lengths = len
  result = v
  while sum_lengths < new_length
    if len + sum_lengths <= new_length
      result = vcat(result,v)
      sum_lengths = sum_lengths+len
    else
      result = vcat(result,v[1:(new_length-sum_lengths)])
      sum_lengths = sum_lengths+new_length-len
    end
    #println("sum_lengths: ",sum_lengths,"  result: ",result)
  end
  result
end


# Works on windows if current working directory (by pwd()) is C:\\Users\\oldmtnbiker\\Dropbox\\evotech\\complexity\\data
function goal_complexity_frequency_dataframe( p::Parameters, gl_init::GoalList, 
    cmplx_df_csvfile::String="../data/consolidate/geno_pheno_raman_df_all_9_13.csv", 
    freq_df_csvfile::String="../data/counts/count_out_4x1_all_ints_10_10.csv" )
  gl_string = map(x->@sprintf("0x%x",x[1]),gl_init) # convert goals from Vector of MyInts to strings
  result_df = DataFrame()
  result_df.goal = gl_init
  cdf = read_dataframe(cmplx_df_csvfile)  # get datafrane with complex and counts
  result_df.complexity =  [cdf[cdf.goal.==gl_string[i],:complex][1] for i = 1:length(gl_string)] # add complexity
  fdf = read_dataframe(freq_df_csvfile)
  result_df.ints8_5 =  [fdf[fdf.goals.==gl_string[i],:ints8_5][1] for i = 1:length(gl_string)] # add complexity
  result_df.ints9_5 =  [fdf[fdf.goals.==gl_string[i],:ints9_5][1] for i = 1:length(gl_string)] # add complexity
  result_df.ints10_5 =  [fdf[fdf.goals.==gl_string[i],:ints10_5][1] for i = 1:length(gl_string)] # add complexity
  result_df.ints11_5 =  [fdf[fdf.goals.==gl_string[i],:ints11_5][1] for i = 1:length(gl_string)] # add complexity
  result_df.ints11_8 =  [fdf[fdf.goals.==gl_string[i],:ints11_8][1] for i = 1:length(gl_string)] # add complexity
  result_df
end
  
# Uses randgoallist to generate a random goal list of length gl_init_length.
# Reads a dataframe that contains complexities for all goals
# Creates a dataframe with goals and complexities, and sorts by complexities
# Returns the list of goals of length gl_final_length of greater complexity.
# Assumes goals have 1 component
# Works on windows if current working directory (by pwd()) is C:\\Users\\oldmtnbiker\\Dropbox\\evotech\\complexity\\data
function filter_goallist_by_complexity( p::Parameters, gl_init::GoalList, gl_final_length::Int64, 
      cmplx_df_csvfile::String="../data/consolidate/geno_pheno_raman_df_all_9_13.csv")
  #gl_init = randgoallist(gl_init_length,p.numinputs, p.numoutputs )
  gl_string = map(x->@sprintf("0x%x",x[1]),gl_init) # convert goals from Vector of MyInts to strings  
  result_df = DataFrame()
  result_df.goal = gl_init
  cdf = read_dataframe( cmplx_df_csvfile ) # get datafrane with :complex field
  result_df.complexity =  [cdf[cdf.goal.==gl_string[i],:complex][1] for i = 1:length(gl_string)] # add complexity
  sort!( result_df, order(:complexity, rev=true))
  result_df[1:gl_final_length,:goal]
end

# Possible figfile:  "11_7/complexities_$(length(gl))_random_goals.png"
function complexity_freq_scatter_plot( p::Parameters, gl::GoalList, inverse_power::Int64, count_field::Symbol=:ints11_8; 
     denom::Int64=2, figfile::String="")
  inverse_root( x, power, denom ) = x^(1/power)/denom
  ngoals = length(gl)
  gcdf = goal_complexity_frequency_dataframe( p, gl )  
  gcdf.crf = map(x->inverse_root(x,inverse_power,denom),gcdf[!,count_field])
  sort!(gcdf, order(:crf, rev=true))  
  rp= Random.shuffle!(collect(1:size(gcdf)[1]))   # random permutation of indices
  root_types = ["$(inverse_power)th","square", "cube", "fourth"]
  root_type = inverse_power in [2,3,4] ? root_types[inverse_power] : root_types[1]
  xlabel = "dot diameter is $(root_type) root of goal frequency"
  title = "Complexities of $ngoals random goals"
  plt = Plots.scatter(rp, gcdf.complexity, markersize=gcdf.crf, label="", grid="none",xticks=[],xshowaxis=false,ylabel="complexity",xlabel=xlabel,markercolor=collect(1:ngoals),title=title)
  if length(figfile) > 0
    savefig(figfile)
  end
  plt
end

#  Bin values greater than or equal to v_min with boundaries [v_min, vmin+1.0/bin_fract_denom, v_min+2.0/bin_fract_denom, ... ]
#  Example:  to bin v into intervals  [2.8, 2.9), [2.9,3.0), [3.0,3.1) ... , using bin_value( v, 2.8, 10.0 )
#  values less than v_min will be mapped to 1.
function bin_value( v::Float64, v_min::Float64, bin_fract_denom::Float64 )  
  result = Int(floor( bin_fract_denom*(v - v_min )))
  result >= 0.0 ? result+1 : 1    # Increase by 1 since Julia uses 1-based indexing
end

# Example:  let values = [0.3, 0.47, 0.57, 0.4, 0.76, 0.67, 0.04, 0.01, 0.13, 0.88]
# The bin_data( values, 10.0) bins values into intervals [0.0,0.1), [0.1,0.2), [0.2,0.3), . . .
# So bin_data( values, 10.0 ) = [2, 1, 0, 1, 2, 1, 1, 1, 1]  # transposed 
function bin_data( values::Vector{Float64}, bin_fract_denom::Float64 )
  v_min = floor(bin_fract_denom*minimum(values))/bin_fract_denom
  ind_max = bin_value(maximum(values),v_min,bin_fract_denom)
  vdict = Dict{Int64,Int64}( i=>0 for i = 1:ind_max ) 
  for v in values
    vdict[bin_value(v, v_min, bin_fract_denom )] += 1
  end
  [ vdict[key] for key in sort(collect(keys(vdict))) ]
end

#function bin_counts( complexity_count_pairs::Vector{Tuple{Float64,Int64}}, bin_fract_denom::Float64 ) 
function bin_counts( values::Vector{Float64}, counts::Vector{Float64}, bin_fract_denom::Float64 ) 
  #values = map(x->x[1],complexity_count_pairs)
  @assert length(values) == length(counts)
  v_min = floor(bin_fract_denom*minimum(values))/bin_fract_denom
  ind_max = bin_value(maximum(values),v_min,bin_fract_denom)
  vdict = Dict{Int64,Int64}( i=>0 for i = 1:ind_max ) 
  for i in 1:length(values)
    vdict[bin_value(values[i], v_min, bin_fract_denom )] += counts[i]
  end
  [ vdict[key] for key in sort(collect(keys(vdict))) ]
end

# Run kolmogorov_complexity() for all goals in goal list gl (in parallel).
# Write dataframe to csvfile is that is given (keyword argument)
function run_kolmogorov_complexity( p::Parameters, funcs::Vector{Func}, gl::GoalList, max_goal_tries::Int64, max_ev_steps::Int64;
      use_mut_evolve::Bool=false, use_lincircuit::Bool=false, decimal_goals::Bool=false, csvfile::String="" ) 
  ngoals = length(gl)
  df = DataFrame()
  df.goal = Goal[]
  df.num_gates = Int64[]
  df.num_active = Int64[]
  df.complexity = Float64[]
  df.tries = Int64[]
  df.avg_robustness = Float64[]
  df.avg_evolvability = Float64[]
  df.num_gates_exc = Int64[]
  result_list = pmap( g->kolmogorov_complexity( p, funcs, g, max_goal_tries, max_ev_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve ), gl )
  #result_list = map( g->kolmogorov_complexity( p, funcs, g, max_goal_tries, max_ev_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve ), gl )
  for r in result_list
    push!(df, r )
  end
  if decimal_goals
    df.goal = map( g->Int(g[1]), df.goal )
  end
  hostname = readchomp(`hostname`)
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# max_goal_tries: ",max_goal_tries)
      println(f,"# max_ev_steps: ",max_ev_steps)
      println(f,"# ngoals: ",ngoals)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  println("# date and time: ",Dates.now())
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )
  println("# funcs: ", Main.CGP.default_funcs(p.numinputs))
  print_parameters(p,comment=true)
  println("# max_goal_tries: ",max_goal_tries)
  println("# max_ev_steps: ",max_ev_steps)
  println("# ngoals: ",ngoals)
  df
end

# Try to find the minimum number of gates to evolve a goal which I call the Kolmogorov complexity.
# Start with p.numinteriors and then deccrease the number of gates until the goal is found.
# Do max_goal_tries evolutions with each number of gates.
# In the unusual case when no circuit is found on this first loop, a dataframe record with num_gates=-1 is returned. 
#   This case should be fixed by hand.
# Then once min number gates is found, do up to num_tries_multiplier*max_goal_tries further evolutions to get
#    max_goal_tries further approximations of Tononi complexity, robustness, and evolvability.
# If a chromosome with a smaller number of active gates is found in these further evolutions,
#    then num_gates is reset and the further evaluations of complexity, robustness, evolvability are restarted
# max_ev_steps is the maximum number of steps while doing mut_evolve()
function kolmogorov_complexity( p::Parameters, funcs::Vector{Func}, g::Goal, max_goal_tries::Int64, max_ev_steps::Int64; use_lincircuit::Bool=false, use_mut_evolve::Bool=false )
  num_tries_multiplier = 3   # used in second outer while loop
  if use_lincircuit && use_mut_evolve
    use_mut_evolve = false   # use_mut_evolve==true is incompatible with use_lincircuit==true
  end
  println("goal: ",g)
  num_gates_exceptions = 0
  #funcs = default_funcs(p.numinputs)
  complexities_list = Float64[]
  #robust_evol_list = Tuple{Float64,Float64}[]
  robust_list = Float64[]
  evol_list = Float64[]
  # if use_lincircuit then num_gates is actually num_registers
  num_gates = p.numinteriors+1   # decrement on the first iteration
  p_current = p
  found_c = use_lincircuit ?  LinCircuit(Vector{UInt16}[[],[],[],[],[],[]],p) : Chromosome(p,[],[],[],0.0,0.0) # dummy circuit to establish scope
  goal_found = true
  while num_gates > 1 && goal_found  # terminates when no goal is found for this value of num_gates
    num_gates -= 1
    #println("while num_gates > 1 && goal_found: num_gates: ",num_gates)
    #println("num_gates: ",num_gates)
    p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
    goal_found = false
    tries = 0
    while !goal_found && tries < max_goal_tries  
      #println("   while !goal_found && tries < max_goal_tries: tries: ",tries,"  num_gates: ",num_gates)
      tries += 1
      #println("inner while tries: ",tries)
      c = use_lincircuit ? rand_lcircuit( p_current, funcs) :  random_chromosome( p_current, funcs )
      if use_mut_evolve
        (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [g], funcs, max_ev_steps, print_steps=false ) 
      else
        (c,step) = neutral_evolution( c, funcs, g, max_ev_steps )
      end
      if step < max_ev_steps
        outputs = output_values( c )
        if sort(outputs) != sort(g)
          println("g: ",g,"  outputs: ",outputs)
        end
        try
          @assert sort(outputs) == sort(g)
        catch(e)
          println("g: ",g,"  outputs: ",outputs)
        end
        goal_found = true
        found_c = deepcopy(c)
        num_gates = number_active(found_c)
        p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
        println("circuit found for goal ",g," with num_gates = ",num_gates,"  output values: ",output_values(found_c))
      end
    end
  end
  if num_gates >= 1 && !goal_found
    println("no goal found for goal ",g," for num_gates: ",num_gates,"  num_gates set to: ", num_gates+1 )
    num_gates += 1  # now set to the minimum number of gates for successful circuit
  end  
  k_failure = use_lincircuit ? length(found_c.circuit_vects[1])==0 : length(found_c.outputs)==0
  if k_failure
    println("found_c not found: ",found_c)
    return (g,-1,-1,-1.0,tries,-1.0,-1.0,num_gates_exceptions)
  end
  if p_current.numinteriors <= 20   # Too time consuming for a large number of gates
    push!(complexities_list, use_lincircuit ? lincomplexity(found_c,funcs) : complexity5(found_c))
  end
  #push!(robust_evol_list, mutate_all( found_c, funcs,robustness_only=true))
  #println("end first loop:  mutate_all ")
  (rlist,elist) = mutate_all( found_c, funcs,robustness_only=true)
  push!(robust_list, rlist )
  push!(evol_list, elist )
  #println("max iter: num_tries_multiplier*max_goal_tries: ",num_tries_multiplier*max_goal_tries)
  p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
  tries = 1
  iter = 0  # put a bound on iterations
  # Iteration bound might be needed because when a circuit is found with fewer active gates than num_gates, tries is reset to 1.
  while iter < num_tries_multiplier*max_goal_tries && tries < max_goal_tries
    print("goal: ",g,"  iter: ",iter,"  tries: ",tries,"  num_gates: ",num_gates,"   ")
    c = use_lincircuit ? rand_lcircuit( p_current, funcs) :  random_chromosome( p_current, funcs )
    println("mut_evolve for goal ",g,"  with numints : ",c.params.numinteriors)
    if use_mut_evolve
      (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [g], funcs, max_ev_steps, print_steps=true ) 
    else
      (c,step) = neutral_evolution( c, funcs, g, max_ev_steps )
    end
    if step < max_ev_steps
      outputs = output_values( c )
      #println("step: ",step,"  outputs: ",outputs,"  goal: ",g)
      @assert sort(outputs) == sort(g)
      if num_gates != number_active(c)  # In this case, tries is reset to 1
        println("num_gates not equal to number active gates for goal: ",g)
        println("num_gates: ",num_gates,"  number_active(c): ",number_active(c))
        #print_build_chromosome(c)
        found_c = c =  remove_inactive( c )
        num_gates = number_active(found_c)
        p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
        num_gates_exceptions += 1
        complexities_list = Float64[]  # restart complexities_list
        #robust_evol_list = Tuple{Float64,Float64}[]  # restart robust_evol_list 
        robust_list = Float64[]  # Restart robust_list
        evol_list = Float64[]    # Restart evol_list
        iter = 0
        tries = 1  # do more tries
      end
      push!(complexities_list, use_lincircuit ? lincomplexity(found_c,funcs) : complexity5(found_c))
      #push!(robust_evol_list, mutate_all( c, funcs,robustness_only=true))
      (rlist,elist) = mutate_all( c, funcs,robustness_only=true)
      push!(robust_list, rlist )
      push!(evol_list, elist )
      tries += 1
    end
    iter += 1
  end
  #println("tries: ",tries)
  avg_complexity = sum(complexities_list)/length(complexities_list)
  #avg_robustness = sum( map(x->x[1],robust_evol_list))/length(robust_evol_list)
  avg_robustness = sum(robust_list)/length(robust_list)
  #avg_evolvability = sum( map(x->x[2],robust_evol_list))/length(robust_evol_list)
  avg_evolvability = sum(evol_list)/length(evol_list)
  (g,num_gates,number_active(found_c),avg_complexity,tries,avg_robustness,avg_evolvability,num_gates_exceptions)
end

function k_complexity_density( p::Parameters, funcs::Vector{Func}, phlist::GoalList; csvfile::String="" )
  kdict = kolmogorov_complexity_dict(p,funcs)
  k_comp_list = Int64[]
  for ph in phlist
    push!(k_comp_list,kdict[ph[1]])
  end
  df = DataFrame( :phlist=>phlist, :num_gates=>k_comp_list )  
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs)
      print_parameters(f,p,comment=true)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function lg_redund_density( p::Parameters, funcs::Vector{Func}, phlist::GoalList; csvfile::String="" )
  rdict = redundancy_dict(p,funcs)
  redund_list = Int64[]
  for ph in phlist
    push!(redund_list,rdict[ph[1]])
  end
  df = DataFrame( :phlist=>phlist, :lg_redund=>map(lg10,redund_list ) )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs)
      print_parameters(f,p,comment=true)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function robustness_density( p::Parameters, funcs::Vector{Func}, phlist::GoalList; csvfile::String="" )
  if p.numinputs == 4 && p.numinteriors == 10 && p.numlevelsback == 5
    rdf = read_dataframe( "../data/12_7_22/ph_evolve_12_7_22G.csv")
  else
    error("only 4x1 10gts 5lb allowed")
  end
  grobust_list = Float64[]
  for i = 1:500
    push!( grobust_list, robustness( random_chromosome( p, funcs), funcs ) )
  end  
  df = DataFrame( :phlist=>phlist, :ph_robust=>rdf.robustness, :geno_robust=>grobust_list )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs)
      print_parameters(f,p,comment=true)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function mutational_density( p::Parameters, funcs::Vector{Func}, ph::Goal, max_tries::Int64, max_steps::Int64 )
  ch = pheno_evolve( p, funcs, ph, max_tries, max_steps )[1]
  alist = map(x->x[1],reduce(vcat,map(x->mutate_all(x),mutate_all(ch,output_circuits=true,output_outputs=false))))
  dlist = map(x->hamming(output_values(ch)[1],x),alist)
  #density(dlist)
  #plot!(title="density_double_mutations_of_circuit_pheno_0xa5a5")
  #savefig("../data/12_31_22/density_double_mutations_of_circuit_pheno_0xa5a5")o
end

# For each phenotypes in goallist and for each numints in the range numinteriors, 
#    computes the mean and standard deviation of the Tononi complexity of numcircuits chromosomes evolved to map to the phenotype.
# Purpose:  try to understand how Tononi complexity scales as numinteriors.
# If normalize==true, then divides the computed Tononi complexity by numints-1 to test for linear scaling.  Does not show linear scaling.
# Tests:  data/7_26_22/
# Assumes 1 output
function tononi_complexity_multiple_params( numinputs::Int64, numlevelsback::Int64, numinteriors::AbstractRange, numcircuits::Int64, goallist::GoalList, max_tries::Int64, max_steps::Int64;
      normalize::Bool=false, csvfile::String="" )
  println("normalize: ",normalize)
  df = DataFrame()
  df.goal = Vector{MyInt}[]
  df.numcircs = Int64[]
  for i in numinteriors
    #println("i: ",i,"  sym: ",Symbol("mean_cmplx_$(i)_ints"))
    insertcols!(df,Symbol("meancmplx$(i)ints")=>Float64[])
  end
  for i in numinteriors
    insertcols!(df,Symbol("stdcmplx$(i)_ints")=>Float64[])
  end
  for ph in goallist
    println("ph: ",ph)
    mean_complexities = zeros(Float64,length(numinteriors))
    std_complexities = zeros(Float64,length(numinteriors))
    i=1
    for numints in numinteriors
      p = Parameters( numinputs, 1, numints, numlevelsback )
      funcs = default_funcs( numinputs )
      circuit_steps_list = pheno_evolve( p, funcs, ph, numcircuits, max_tries, max_steps )
      circuits_list = map(x->x[1],circuit_steps_list)
      steps_list = map(x->x[2],circuit_steps_list)
      complexity_list = map(c->complexity5(c), circuits_list )
      println("numints: ",numints, "  complexity_list: ",complexity_list)
      if normalize
        complexity_list ./= (numints-1)
      end
      println("numints: ",numints, "  complexity_list: ",complexity_list)
      mean_complexities[i] = mean(complexity_list)
      std_complexities[i] = std(complexity_list)
      i += 1
    end
    df_row = vcat([ ph, numcircuits], mean_complexities, std_complexities )
    println("size(df): ",size(df),"  length(df_row): ",length(df_row))
    push!( df, df_row )
  end
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(numinputs))
      println(f,"# numinputs: ",numinputs)
      println(f,"# numlevelsback: ",numlevelsback)
      println(f,"# numinteriors: ",numinteriors)
      println(f,"# normalize: ",normalize)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function redund_vs_k_complexity_plot( p::Parameters, funcs::Vector{Func}, phlist::GoalList, title::String="K complexity vs log redundancy ";
      csvfile::String="", plotfile::String="" )
  lg10(x::Number) = x > 0 ? log10(x) : 0
  k_dict = kolmogorov_complexity_dict( p )
  r_dict = redundancy_dict( p, csvfile )
  ttitle = string(title,"$(p.numinputs)x1 $(p.numinteriors)gts $(p.numlevelsback)lb $(length(funcs))funcs")
  println("ttitle: ",ttitle)
  goals = map(x->x[1],phlist)
  kcomp = map(x->k_dict[x],goals)
  log_redund = map(x->lg10(r_dict[x]),goals)
  pdf = DataFrame( :goal=>map(x->@sprintf("0x%04x",x),goals), :kcomp=>kcomp, :log_redund=>log_redund )
  gr()
  plt = scatter(pdf.kcomp,pdf.log_redund,labels=:none,xlabel="Kolmogorov complexity",ylabel="log redundancy",title=ttitle)
  if plotfile != ""
    savefig(plotfile)
  end
  #pdf
end

function kolmogorov_complexity_dict( p::Parameters, funcs::Vector{Func}=default_funcs(p), csvfile::String="" )
  if p.numinputs == 3
    if length(funcs) == 4
      k_csvfile = "../data/counts/k_complexity_3x1_4funcs7_11_22E.csv"
    elseif length(funcs) == 5
      k_csvfile = "../data/counts/k_complexity_3x1_5funcs7_11_22F.csv"
    end
  elseif p.numinputs == 4
    if length(funcs) == 5
      #k_csvfile = "../data/counts/k_complexity_all4x1phenos.csv"
      k_csvfile = "../data/counts/k_complexity8_9_22FGGF.csv"
    else
      k_csvfile = "../data/counts/k_complexity8_9_22XTYZZYTZ.csv"   
    end
  else
    println("Illegal parameters or funcs in call to kolmogorov_complexity_dict() ")
    return nothing
  end
  println("k_csvfile: ",k_csvfile)
  df = read_dataframe( k_csvfile )
  dict = Dict{ MyInt, Int64 }()
  pheno_name = Symbol(names(df)[1])
  #println("size(df): ",size(df))
  for i = 1:size(df)[1]
    dict[string_to_MyInt(df[i,pheno_name])] = df.num_gates[i]
  end
  dict
end

function redundancy_dict( p::Parameters, funcs::Vector{Func}=default_funcs(p), csvfile::String="" )
  if length(csvfile) == 0
    if p.numinputs == 3
      if p.numinteriors == 8 && p.numlevelsback == 4
        if length(funcs) == 4
          csvfile = "../data/counts/count_outputs_ch_4funcs_3inputs_8gates_4lb_W.csv"
        elseif length(funcs) == 5
          csvfile = "../data/counts/count_outputs_ch_5funcs_3inputs_8gates_4lb_V.csv"
        else
          println("no csv file in function redundancy_dict()")
          return nothing
        end
      elseif p.numinteriors == 7 && p.numlevelsback == 4
        if length(funcs) == 4
          csvfile = "../data/counts/count_outputs_ch_4funcs_3inputs_7gates_4lb_cmplxC.csv"
        else
          println("no csv file in function redundancy_dict()")
          return nothing
        end
      elseif p.numinteriors == 8 && p.numlevelsback == 5
        if length(funcs) == 4
          csvfile = "../data/counts/count_outputs_3x1_8gts5lb_4funcs.csv"
        else
          println("no csv file in function redundancy_dict()")
          return nothing
        end
      elseif p.numinteriors == 10 && p.numlevelsback == 5
        if length(funcs) == 5
          csvfile = "../data/counts/count_outputs_ch_5funcs_3inputs_10gates_5lb_Y.csv"
        else
          println("no csv file in function redundancy_dict()")
          return nothing
        end
      else
        println("no csv file in function redundancy_dict()")
        return nothing
      end
    elseif p.numinputs == 4
      if length(funcs) == 5 
        if p.numinteriors == 10 && p.numlevelsback == 5
          csvfile = "../data/counts/count_outputs_ch_5funcs_4inputs_10gates_5lb_EG.csv"
        elseif p.numinteriors == 12 && p.numlevelsback == 6
          csvfile = "../data/counts/count_outputs_ch_5funcs_4inputs_12gates_6lb_H.csv"
        elseif p.numinteriors == 11 && p.numlevelsback == 8
          csvfile = "../data/counts/count_out_4x1_all_ints_11_8.csv"
        else
          println("only 10gts 5lb and 12gts 6 lb are supported at this time")
          return nothing
        end
      else
        println("Only 5 funcs is supported for 4 inputs at this time")
        return nothing
      end
    else
      println("only 3 and 4 inputs are supported at this time")
      return nothing
    end
  end
  println("csvfile: ",csvfile)
  df = read_dataframe( csvfile )
  dict = Dict{ MyInt, Int64 }()
  pheno_name = Symbol(names(df)[1])
  counts_name = Symbol(names(df)[2])
  #println("pheno_name: ",pheno_name,"  counts_name: ",counts_name)
  println("size(df): ",size(df))
  for i = 1:size(df)[1]
    dict[string_to_MyInt(df[i,pheno_name])] = df[i,counts_name]
  end
  dict
end

function randunique( lst::Vector, numelements::Int64, numtries::Int64 )
  if length(lst) <= numelements
    return lst
  else
    return_list = type[]
    successes = 0
    iter = 1
    while iter < numtries && successes < numelements
      elt = rand(lst)
      #println("iter: ",iter,"  elt: ",elt, "  return_list: ",return_list)
      iter += 1
      if !(elt in return_list)
        push!(return_list,elt)
        successes += 1
      end
    end
    if successes == numelements
      return return_list
    else
      error("function randunique() failed.  You can try increasing the numtries parameter.")
    end
  end
end

# Given a dataframe whose first column is a column of phenotypes (goals), generate the dataframe with this column replaced by the ones complement of the phenotypes.
# Assumes that for the other columns the value of a row corresponding to a phenotype is the same as the value corresponding to the ones complement of that phenotype
# Parameters p is only used to call default_funcs(p), so only p.numinputs is used.
function negate_dataframe( df::DataFrame, p::Parameters; combine::Bool=false )
  default_funcs(p)  # set global variable Ones so that CGP.Not will work
  m_one = MyInt(1)
  if typeof( df[1,1] ) <: AbstractString
    goals = map( g->eval(Meta.parse( g ))[1], df[:,1] )
  else
    goals = Vector{MyInt}(df[:,1])
  end
  neg_goals = reverse( map( g->[CGP.Not(g)], goals ) )
  ndf = DataFrame( Symbol(names(df)[1])=>neg_goals )
  for i = 2:size(df)[2]
    insertcols!(ndf,Symbol(names(df)[i])=>reverse(df[:,i]))
  end
  if !combine
    return ndf
  else
    ndf.goal = map(g->String15(@sprintf("%s",g)),ndf.goal)
    return vcat( df, ndf )
  end
end

# Prints the phenotypes where the K complexity of every phenotype is not equal to the K complexity of its ones complement
# If the second argument is true, then sets entries where there is a disagreement to the smaller value
# Kcomplx must be a vector indexed over all phenotypes for some number of inputs
function K_complexity_negation_check( Kcmplx::Vector{Int64}, correct::Bool=false )
  count_errors = 0
  m_one = MyInt(1)
  for i = MyInt(0):MyInt(div( length(Kcmplx), 2 ) -1 )
    #println("i+1: ",i+m_one,"  CGP.Not(i)+m_one: ",CGP.Not(i)+m_one,"   ")
    #println("i+1: ",i+1,"  CGP.Not(i)+1: ",CGP.Not(i)+1,"   ")
    #if Kcmplx[i+m_one] != Kcmplx[ CGP.Not(i)+m_one ]
    if Kcmplx[i+1] != Kcmplx[ CGP.Not(i)+1 ]
      count_errors += 1
      #@printf("0x%04x   0x%04x:  ", i+m_one, CGP.Not(i) )
      @printf("0x%04x   0x%04x:  ", i, CGP.Not(i) ) #println(Kcmplx[i+m_one],"  ",Kcmplx[ CGP.Not(i)+m_one ])
      println(Kcmplx[i+1],"  ",Kcmplx[ CGP.Not(i)+1 ])
      if correct
        #knew = min(Kcmplx[i+m_one],Kcmplx[ CGP.Not(i)+m_one ])
        knew = min(Kcmplx[i+1],Kcmplx[ CGP.Not(i)+1 ])
        #Kcmplx[i+m_one] = knew
        Kcmplx[i+1] = knew
        #Kcmplx[ CGP.Not(i)+m_one ] = knew
        Kcmplx[ CGP.Not(i)+1 ] = knew
      end
    end
  end
  count_errors
end

# Prints the phenotypes where the K complexity of every phenotype is not equal to the K complexity of its ones complement
# If the second argument is true, then sets entries where there is a disagreement to the smaller value
# Kcmplx1 and Kcmplx2 must be a vectors indexed over the same collection of phenotypes for some number of inputs
function Kcomplexity_accuracy_check( Kcmplx1, Kcmplx2, correct::Bool=false )
  @assert length(Kcmplx1)==length(Kcmplx2)
  max_error = 4
  count_errors_plus = zeros(Int64,max_error)
  count_errors_minus = zeros(Int64,max_error)
  m_one = MyInt(1)
  for i = MyInt(0):MyInt(length(Kcmplx1)-1)
    if Kcmplx1[i+1] != Kcmplx2[i+1]
      @printf("0x%04x", i )
      println("  ",Kcmplx1[i+1],"  ",Kcmplx2[i+1])
      if Kcmplx1[i+1] > Kcmplx2[i+1]
        count_errors_plus[ Kcmplx1[i+1] - Kcmplx2[i+1] ] += 1
      elseif Kcmplx2[i+1] > Kcmplx1[i+1]
        count_errors_minus[ Kcmplx2[i+1] - Kcmplx1[i+1] ] += 1
      end
      if correct
        k1 = Kcmplx1[i+1]
        k2 = Kcmplx2[i+1]
        if k1 < k2
          Kcmplx2[i+1] = k1
        elseif k2 < k1
          Kcmplx1[i+1] = k2
        end
      end
    end
  end
  ( sum(count_errors_plus)+sum(count_errors_minus), count_errors_plus, count_errors_minus )
end

function run_K_complexity_mutate_all( p::Parameters, funcs::Vector{Func}, circuit_ints_df::DataFrame=DataFrame(); 
     use_lincircuit=false, use_mut_evolve::Bool=false, print_steps::Bool=false, csvfile::String="" )
  run_k_complexity_mutate_all( p, funcs, Goal[], 0, 0, 0, circuit_ints_df=circuit_ints_df, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps )
end

function run_k_complexity_mutate_all( p::Parameters, funcs::Vector{Func}, phlist::GoalList, numcircuits::Int64, max_tries::Int64, max_steps::Int64;
    circuit_ints_df::DataFrame=DataFrame(), use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false, csvfile::String="" )
  k_dict = kolmogorov_complexity_dict( p, funcs )
  df = DataFrame( :goal=>Goal[], :rebased_vect=>RebasedVector[] )
  circuit_ints_list = Vector{Int128}[]   # establish scope
  sampling = false
  println("size(circuit_ints_df)[1]: ",size(circuit_ints_df)[1])
  if size(circuit_ints_df)[1] > 0
    println("length(circuit_ints_df.goals): ",length(circuit_ints_df.goal))
    circuit_ints_df.goal = map(x->[eval(Meta.parse(x))],circuit_ints_df.goals)
    println("length(circuit_ints_df.goal): ",length(circuit_ints_df.goal))
    circuit_ints_df.circuit_ints_list = pmap(x->eval(Meta.parse(x)),circuit_ints_df.circuits_list)
    println("length(circuit_ints_df.circuits_list): ",length(circuit_ints_df.circuits_list))
    phlist = circuit_ints_df.goal
    circuit_ints_list = circuit_ints_df.circuit_ints_list
    sampling = true
  end
  println("length(circuit_ints_list): ",length(circuit_ints_list))
  k_comp_rebased = RebasedVector[]
  #cilist = length(circuit_ints_list)>0 ? 
  result_list = pmap( i->k_complexity_mutate_all( p, funcs, phlist[i], numcircuits, max_tries, max_steps, k_dict, 
      circuit_ints_list=length(circuit_ints_list)>0 ? circuit_ints_list[i] : Int128[],
      use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), 1:length(phlist) )
  #result_list = map( i->k_complexity_mutate_all( p, funcs, phlist[i], numcircuits, max_tries, max_steps, k_dict, 
  #    circuit_ints_list=length(circuit_ints_list)>0 ? circuit_ints_list[i] : Int128[],
  #    use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), 1:length(phlist) )
  for res in result_list
    push!(k_comp_rebased,RebasedVector(res[2],res[3]))
    push!(df,(res[1],RebasedVector(res[2],res[3])))
  end
  sdf = comp_summary_dataframe( df )
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      println(f,"# ",(use_lincircuit ? "LGP" : "CGP"))
      println(f,"# ",(sampling ? "sampling" : "evolution"), " evolvability")
      print_parameters(f,p,comment=true)
      println(f,"# numcircuits: ",numcircuits)
      println(f,"# ngoals: ",length(phlist))
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, sdf, append=true, writeheader=true )
    end
  end
  sdf
end

# Find the K complexity of all genotypes generated by mutate_all() applied to numcircuits genotypes that map to phenotype ph
#    and accumulates their counts in k_complexity_counts.
# The default case is that the genotypes that map to ph are generated by calling pheno_evolve().
# The alternate case is that circuit_ints_list supplies a list of circuit ints of circuits that map to ph.
function k_complexity_mutate_all( p::Parameters, funcs::Vector{Func}, ph::Goal, numcircuits::Int64, max_tries::Int64, max_steps::Int64, k_dict::Dict{MyInt,Int64}; 
     circuit_ints_list::Vector{Int128}=Int128[], use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false )
  default_funcs(p)  # Set global variable Ones on this process
  k_complexity_counts = zeros(Int64,12)  # 12 is an upper bound for possible k_complexities
  k_comp_ph = k_dict[ph[1]]
  if length(circuit_ints_list) == 0
    circuit_steps_list = pheno_evolve( p, funcs, ph::Goal, numcircuits, max_tries::Int64, max_steps::Int64,
        use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps )
    circuit_list = map( x->x[1], circuit_steps_list )
  else
    circuit_list = map( ci->int_to_chromosome(ci,p,funcs),circuit_ints_list)
    println("ph: ",ph,"  length(circuit_list): ",length(circuit_list))
  end
  for c in circuit_list
    outputs_list = mutate_all( c, funcs, output_outputs=true )
    k_comp_list = map( x->k_dict[x[1]], outputs_list )
    for k_comp in k_comp_list
      k_complexity_counts[k_comp] += 1
    end
  end
  (ph,k_comp_ph,k_complexity_counts)
end
    
function comp_summary_dataframe( comp_rebased::DataFrame )
  sdf = DataFrame( :kcomp=>Int64[] ) 
  for i = -8:8
    insertcols!(sdf,Symbol(string(i))=>Int64[])
  end
  println("size(sdf): ",size(sdf))
  summary = [ RebasedVector(u,zeros(Int64,17)) for u = 1:8 ]
  for i = 1:size(comp_rebased)[1]
    j = comp_rebased.rebased_vect[i].center
    for k = -8:8
      summary[j][k] = summary[j][k]+comp_rebased.rebased_vect[i][k] 
      #println("i: ",i,"  j: ",j,"  k: ",k,"  ksummary[j][k]: ",ksummary[j][k])
    end
  end
  for u = 1:8
    #println("row length: ",length(vcat([u],summary[u].vect)))
    push!( sdf, vcat([u],[summary[u][k] for k = -8:8]) )
  end
  sdf
end

# Convert the string column circuit_ints_df.goals to the GoalList column circuit_ints_df.goal.
# And convert the string column circuit_ints_df.circuits_list to the Vector{Int128} column circuit_ints_df.circuit_ints_list.
function convert_circuit_ints_df( circuit_ints_df::DataFrame )
  circuit_ints_df.goal = map(x->[eval(Meta.parse(x))],circuit_ints_df.goals) 
  circuit_ints_df.circuit_ints_list = pmap(x->eval(Meta.parse(x)),circuit_ints_df.circuits_list)
end

function Base.getindex( rb::RebasedVector, i::Integer )
  try  # if bounds error return 0
    rb.vect[i+rb.center]
  catch
    0
  end
end
    
function Base.setindex!( rb::RebasedVector, v::Vector{Int64}, i::Vector{Integer} )
  rb.vect[i+rb.center] = v
end
    
function Base.setindex!( rb::RebasedVector, v::Int64, i::Integer )
  try  # if bounds error do nothing
    rb.vect[i+rb.center] = v
  catch
  end
end

# Computes the mutation evolvability matrix of the results of mutate_all() to all of the circuits given in the circuits_list
#   column of the cnt_csvfile dataframe which is determined at the start of the function.
# E[s+1,d+1]  is the number of mutations from phenotype s that give phenotype d.
# The arguments source and dest specify K complexity values for rows and columns respectively.
# Then phenotypes selected by the ranges source and dest to produce a submatrix of E corresponding to rows in source
#   and columns in dest.
# Returns a dataframe corresponding to this submatrix.
# Objective:  Show that mutations from low complexity phenotypes to high complexity phenotypes decrease exponentially with the difference in complexity.
function K_complexity_mutation( p::Parameters, funcs::Vector{Func}, source::UnitRange, dest::UnitRange; csvfile::String="" ) 
  # Assumes MyInt == UInt16
  if p.numinputs == 2
    phlist = map(x->[x],0x0000:0x000f)
  elseif p.numinputs == 3
    phlist = map(x->[x],0x0000:0x00ff)
    if length(funcs) == 4
      if p.numinteriors == 7
        cnt_csvfile = "../data/counts/count_outputs_ch_4funcs_3inputs_7gates_4lb_cmplxC.csv"
      elseif p.numinteriors == 8
        cnt_csvfile = "../data/counts/count_outputs_ch_4funcs_3inputs_8gates_4lb_W.csv"
      else
        error("no csv file")
      end
    elseif length(funcs) == 5
      cnt_csvfile = "../data/couunts/count_outputs_ch_5funcs_3inputs_8gates_4lb_V.csv"
    end
  elseif p.numinputs == 4
    phlist = map(x->[x],0x0000:0xffff)
    error("no csv file")
  else
    error("no csv file")
  end
  println("cnt_csvfile: ",cnt_csvfile)
  wdf = read_dataframe( cnt_csvfile )
  kdict = kolmogorov_complexity_dict(p,funcs)
  wsdf = evolvable_pheno_df( p, funcs, phlist, wdf.circuits_list )
  insertcols!(wsdf,:k_complexity=>map(x->kdict[x[1]],phlist))
  println("size(wsdf): ",size(wsdf))
  E = pheno_vects_to_evolvable_matrix( wsdf.pheno_vects )
  source_bv = BitVector(map( i->wsdf.k_complexity[i] in collect(source),1:256))
  source_str = wsdf[source_bv,:pheno_list]
  source_list = map(x->eval(Meta.parse(x)),wsdf[source_bv,:pheno_list])
  dest_bv = BitVector(map( i->wsdf.k_complexity[i] in collect(dest),1:256))
  dest_str = wsdf[dest_bv,:pheno_list]
  dest_list = map(x->eval(Meta.parse(x)),wsdf[dest_bv,:pheno_list])
  sddf = submatrix_to_dataframe( p, funcs, E, wsdf, source_list, dest_list )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs )
      println(f,"# source: ",source)
      println(f,"# dest: ",dest)
      println(f,"# matrix nonzeros: ",matrix_dataframe_count_nonzeros(sddf))
      CSV.write( f, sddf, append=true, writeheader=true )
    end
  end
end

# Counts the number of nonzeros in the matrix corresponding to the dataframe df.
# The first column of df should be the phenotypes, and the remaining should be matrix entries
function matrix_dataframe_count_nonzeros( df::DataFrame )
  reduce(+,map( i->sum(df[:,i]), 2:size(df)[2] ))
end
using DataFrames, Combinatorics, CSV

# maxreps is the number of iterations
# maxtries is how often to keep trying when an iteration fails.
# maxsteps is the maximum number of steps in mut_evolve.  Should be large, like =50000 for 3 inputs
function neighbor_complexity( g::Goal, p::Parameters, maxreps::Int64, maxsteps::Int64, maxtries::Int64=2*maxreps )
  funcs = default_funcs(p.numinputs)      
  chrome_complexity_sum = 0.0
  neighbor_complexity_sum = 0.0
  c_count = 0
  n_count = 0
  i = 0
  nrepeats = 0
  while i < maxtries && nrepeats < maxreps
    i += 1
    c = random_chromosome(p,funcs)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
      mut_evolve( c, [g], funcs, maxsteps )
    if step == maxsteps
      println("mut evolve failed for goal: ",g)
      continue
    end
    #print("c: ")
    #print_build_chromosome(c)
    #println("complexity: ",complexity5(c))
    chrome_complexity_sum += complexity5(c)
    c_count += 1
    c_output = output_values(c)
    #println("goal: ",g,"  i: ",i,"  nrepeats: ",nrepeats,"  c_output: ",c_output)
    @assert sort(c_output) == sort(g)   
    output_chromes = mutate_all( c, funcs, output_circuits=true, output_outputs=false )
    for out_c in output_chromes
      #print("out_c: ")
      #print_build_chromosome(out_c[2])
      #println("complexity: ",complexity5(out_c))
      neighbor_complexity_sum += complexity5(out_c)
      n_count += 1
    end
    #print("i: ",i,"  c_count: ",c_count,"  c_complex_sum: ",chrome_complexity_sum)
    #println("  n_count: ",n_count,"  n_complex_sum: ",neighbor_complexity_sum)
    nrepeats += 1
  end
  (g, chrome_complexity_sum/c_count, neighbor_complexity_sum/n_count)
end

function run_neighbor_complexity( ngoals::Int64, p::Parameters, maxreps::Int64, maxsteps::Int64, maxtries::Int64=2*maxreps; csvfile::String="")
  goallist = randgoallist( ngoals, p.numinputs, p.numoutputs )
  run_neighbor_complexity(goallist, p, maxreps,maxsteps,maxtries, csvfile=csvfile )
end
  
function run_neighbor_complexity( goallist::Vector{Vector{MyInt}}, p::Parameters, maxreps::Int64, maxsteps::Int64, maxtries::Int64=2*maxreps; csvfile::String="")
  df = DataFrame()
  df.goal = Goal[]
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numinteriors = Int64[]
  df.numlevelsback = Int64[]
  df.circuit_complexity = Float64[]
  df.neighbor_complexity = Float64[]
  tuple_list = pmap( g->neighbor_complexity( g, p, maxreps, maxsteps, maxtries ), goallist )
  #tuple_list = map( g->neighbor_complexity( g, p, maxreps, maxsteps, maxtries ), goallist )
  for tpl in tuple_list
    row_tuple = (tpl[1], p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, tpl[2], tpl[3])
    push!(df, row_tuple) 
  end
  df
  hostname = readchomp(`hostname`)
  open(csvfile,"w") do f
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    #println(f,"# run time in minutes: ",ttime/60)
    println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))    
    print_parameters(f,p,comment=true)
    println(f,"# maxreps: ",maxreps)
    println(f,"# maxsteps: ",maxsteps)
    println(f,"# maxtries: ",maxtries)
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end         
    
# Test whether any phenotypes that I categorized as K complexity Kmin can be produced with a LGP circuit with 5 instructdions
#= Setup
Kmin = 5
kdf = read_dataframe("../data/8_9_22/k_complexity_LGPsummary.csv");
p = Parameters(3,1,Kmin-1,2);funcs=default_funcs(p)[1:4]
f5 = map(x->x-1,findall(x->x>=Kmin,kdf.K_complexity));
bs = BitSet(f5);
@time testKall( p, funcs, bs )  # 360 seconds on Mac with 8 processes
UInt16[] 
=#
function testKall( p::Parameters, funcs::Vector{Func}, bs::BitSet )
  function tf( rng::UnitRange{Int64} )
    result_list = MyInt[]
    for i in rng
      ov = output_values(circuit_int_to_circuit( i, p, funcs ))[1]
      if ov in bs
        push!(result_list, ov )
      end
    end
    result_list
  end
  niters = div( 200^p.numinteriors, 25 )  # Exact result, no remainder
  ranges = [ (i-1)*niters:i*niters for i = 1:25 ]
  results = pmap( i->tf(ranges[i]), 1:25 )
  reduce( vcat, results )
end

function run_K_complexity_mutate_all( p::Parameters, funcs::Vector{Func}, circuit_ints_df::DataFrame=DataFrame(); 
     use_lincircuit=false, use_mut_evolve::Bool=false, print_steps::Bool=false, csvfile::String="" )
  run_redundancy_mutate_all( p, funcs, Goal[], 0, 0, 0, circuit_ints_df=circuit_ints_df, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps )
end

function run_redundancy_mutate_all( p::Parameters, funcs::Vector{Func}, phlist::GoalList, numcircuits::Int64, max_tries::Int64, max_steps::Int64;
    circuit_ints_df::DataFrame=DataFrame(), use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false, csvfile::String="" )
  r_dict = redundancy_dict( p, funcs )
  df = DataFrame( :goal=>Goal[], :rebased_vect=>RebasedVector[] )
  circuit_ints_list = Vector{Int128}[]   # establish scope
  sampling = false
  #println("size(circuit_ints_df)[1]: ",size(circuit_ints_df)[1])
  if size(circuit_ints_df)[1] > 0
    println("length(circuit_ints_df.goals): ",length(circuit_ints_df.goal))
    circuit_ints_df.goal = map(x->[eval(Meta.parse(x))],circuit_ints_df.goals)
    println("length(circuit_ints_df.goal): ",length(circuit_ints_df.goal))
    circuit_ints_df.circuit_ints_list = pmap(x->eval(Meta.parse(x)),circuit_ints_df.circuits_list)
    println("length(circuit_ints_df.circuits_list): ",length(circuit_ints_df.circuits_list))
    phlist = circuit_ints_df.goal
    circuit_ints_list = circuit_ints_df.circuit_ints_list
    sampling = true
  end
  #println("length(circuit_ints_list): ",length(circuit_ints_list))
  r_comp_rebased = RebasedVector[]
  result_list = pmap( i->redundancy_mutate_all( p, funcs, phlist[i], numcircuits, max_tries, max_steps, r_dict, 
      circuit_ints_list=length(circuit_ints_list)>0 ? circuit_ints_list[i] : Int128[],
      use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), 1:length(phlist) )
  #result_list = map( i->redundancy_mutate_all( p, funcs, phlist[i], numcircuits, max_tries, max_steps, r_dict, 
  #    circuit_ints_list=length(circuit_ints_list)>0 ? circuit_ints_list[i] : Int128[],
  #    use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), 1:length(phlist) )
  for res in result_list
    push!(r_comp_rebased,RebasedVector(res[2],res[3]))
    push!(df,(res[1],RebasedVector(res[2],res[3])))
  end
  sdf = rcomp_summary_dataframe( df )
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      println(f,"# ",(use_lincircuit ? "LGP" : "CGP"))
      println(f,"# ",(sampling ? "sampling" : "evolution"), " evolvability")
      print_parameters(f,p,comment=true)
      println(f,"# numcircuits: ",numcircuits)
      println(f,"# ngoals: ",length(phlist))
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, sdf, append=true, writeheader=true )
    end
  end
  sdf
end

# Find the K complexity of all genotypes generated by mutate_all() applied to numcircuits genotypes that map to phenotype ph
#    and accumulates their counts in redundancy_counts.
# The default case is that the genotypes that map to ph are generated by calling pheno_evolve().
# The alternate case is that circuit_ints_list supplies a list of circuit ints of circuits that map to ph.
function redundancy_mutate_all( p::Parameters, funcs::Vector{Func}, ph::Goal, numcircuits::Int64, max_tries::Int64, max_steps::Int64, r_dict::Dict{MyInt,Int64}; 
     circuit_ints_list::Vector{Int128}=Int128[], use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false )
  default_funcs(p)  # Set global variable Ones on this process
  redundancy_counts = zeros(Int64,14)  # 14 is an upper bound for possible log redundancies
  r_comp_ph = Int(round(lg10(r_dict[ph[1]])))
  if length(circuit_ints_list) == 0
    circuit_steps_list = pheno_evolve( p, funcs, ph::Goal, numcircuits, max_tries::Int64, max_steps::Int64,
        use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps )
    circuit_list = map( x->x[1], circuit_steps_list )
  else
    circuit_list = map( ci->int_to_chromosome(ci,p,funcs),circuit_ints_list)
    println("ph: ",ph,"  length(circuit_list): ",length(circuit_list))
  end
  for c in circuit_list
    outputs_list = mutate_all( c, funcs, output_outputs=true )
    r_comp_list = map( x->Int(round(lg10(r_dict[x[1]]))), outputs_list )
    r_comp_list = map( r_comp->( iszero(r_comp) ? 1 : r_comp), r_comp_list )
    for r_comp in r_comp_list
      if iszero(r_comp)
        r_comp = 1
        println("rcomp zero r_comp_list: ",r_comp_list)
      end
      redundancy_counts[r_comp] += 1
    end
  end
  (ph,r_comp_ph,redundancy_counts)
end
    
function rcomp_summary_dataframe( r_comp_rebased::DataFrame )
  sdf = DataFrame( :rcomp=>Int64[] ) 
  for i = -8:8
    insertcols!(sdf,Symbol(string(i))=>Int64[])
  end
  println("size(sdf): ",size(sdf))
  rsummary = [ RebasedVector(u,zeros(Int64,17)) for u = 1:8 ]
  for i = 1:size(r_comp_rebased)[1]
    j = r_comp_rebased.rebased_vect[i].center
    for k = -8:8
      try
        rsummary[j][k] = rsummary[j][k]+r_comp_rebased.rebased_vect[i][k] 
      catch
        println("rebased vector error: i: ",i,"  k: ",k,"  j: ",j)
      end
      #println("i: ",i,"  j: ",j,"  k: ",k,"  rsummary[j][k]: ",rsummary[j][k])
    end
  end
  for u = 1:8
    #println("row length: ",length(vcat([u],rsummary[u].vect)))
    push!( sdf, vcat([u],[rsummary[u][k] for k = -8:8]) )
  end
  sdf
end
