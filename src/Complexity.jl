# Includes explore_complexity()
using DataFrames, Statistics, Dates, CSV, Distributed, Printf, Random
export evolve_complexity, run_complexity_list, complexity_list
export run_explore_complexity, explore_complexity
export average_dataframes, extend_list_by_dups, goal_complexity_frequency_dataframe
export filter_goallist_by_complexity, complexity_freq_scatter_plot
export bin_value, bin_data, bin_counts
export kolmogorov_complexity, run_kolmogorov_complexity

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
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
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
  result_dfs = pmap(x->run_explore_complexity( nruns, p, goallist, num_circuits, max_ev_steps,
        insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob), collect(1:nreps))
  #result_dfs = map(x->run_explore_complexity( nruns, p, goallist, num_circuits, max_ev_steps, 
  #      insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob), collect(1:nreps))
  println("length(result_dfs): ",length(result_dfs))
  df = average_dataframes( result_dfs )
  insertcols!(df, 1, :generation=>collect(1:size(df)[1]))
  cumm_unique_goals = zeros(Float64,size(df)[1])
  cumm_unique_goals[1] = df.unique_goals[1]
  for i = 2:size(df)[1]
    cumm_unique_goals[i] = df.unique_goals[i] + cumm_unique_goals[i-1]
  end
  insertcols!(df, size(df)[2]+1, :cumm_unique_goals=>cumm_unique_goals )
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
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
  df
end

# Run explore_complexity() nruns times with the discovered circuit_list from one run becoming the starting circuit_list for the next.
# Also, discoverd goals are removed from goallist on subsequent runs.
function run_explore_complexity( nruns::Int64, p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64;
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  df = DataFrame()
  done = false
  explore_complexity_tries = 0
  while !done
    explore_complexity_tries += 1
    #println("while !done loop.")
    df = DataFrame()
    df.numgates = Float64[]
    df.levelsback = Float64[]
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
  
# This is the explore_complexity() function that runs one generation starting with ciruits in circuit_list.
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
function run_kolmogorov_complexity( p::Parameters, gl::GoalList, max_goal_tries::Int64, max_ev_steps::Int64;
      csvfile::String="" ) 
  ngoals = length(gl)
  df = DataFrame()
  df.goal = Goal[]
  df.num_gates = Int64[]
  df.num_active_gates = Int64[]
  df.complexity = Float64[]
  df.tries = Int64[]
  df.avg_robustness = Float64[]
  df.avg_evolvability = Float64[]
  df.num_gates_exc = Int64[]
  result_list = pmap( g->kolmogorov_complexity( p, g, max_goal_tries, max_ev_steps ), gl )
  #result_list = map( g->kolmogorov_complexity( p, g, max_goal_tries, max_ev_steps ), gl )
  for r in result_list
    push!(df, r )
  end
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
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
# Then once min number gates is found, do up to num_tries_multiplier*max_goal_tries further evolutions to get
#    max_goal_tries further approximations of Tononi complexity, robustness, and evolvability.
# If a chromosome with a smaller number of active gates is found in these further evolutions,
#    then num_gates is reset and the further evaluations of complexity, robustness, evolvability are restarted
# max_ev_steps is the maximum number of steps while doing mut_evolve()
function kolmogorov_complexity( p::Parameters, g::Goal, max_goal_tries::Int64, max_ev_steps::Int64 )
  num_tries_multiplier = 3   # used in second outer while loop
  println("goal: ",g)
  num_gates_exceptions = 0
  funcs = default_funcs(p.numinputs)
  complexities_list = Float64[]
  robust_evol_list = Tuple{Float64,Float64}[]
  num_gates = p.numinteriors+1   # decrement on the first iteration
  p_current = p
  found_c = Chromosome(p,[],[],[],0.0,0.0)  # dummy chromosome to establish scope
  goal_found = true
  while num_gates > 1 && goal_found  # terminates when no goal is found for this value of num_gates
    num_gates -= 1
    #println("num_gates: ",num_gates)
    p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
    goal_found = false
    tries = 0
    while !goal_found && tries < max_goal_tries  
      tries += 1
      #println("inner while tries: ",tries)
      c = random_chromosome( p_current, funcs )
      (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
          mut_evolve( c, [g], funcs, max_ev_steps, print_steps=false ) 
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
        num_gates = number_active_gates(found_c)
        p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
        println("ciruit found for goal ",g," with num_gates = ",num_gates)
      end
    end
  end
  if num_gates >= 1 && !goal_found
    println("no goal found for goal ",g," for num_gates: ",num_gates)
    num_gates += 1  # now set to the minimum number of gates for successful circuit
  end  
  push!(complexities_list, complexity5(found_c))
  push!(robust_evol_list, mutate_all( found_c, funcs,robustness_only=true))
  p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
  tries = 1
  iter = 0  # put a bound on iterations
  while iter < num_tries_multiplier*max_goal_tries && tries < max_goal_tries
    c = random_chromosome( p_current, funcs )
    #println("mut_evolve for goal ",g,"  with numints : ",c.params.numinteriors)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
        mut_evolve( c, [g], funcs, max_ev_steps, print_steps=true ) 
    if step < max_ev_steps
      outputs = output_values( c )
      #println("outputs: ",outputs,"  goal: ",g)
      @assert sort(outputs) == sort(g)
      if num_gates != number_active_gates(c)
        println("num_gates: ",num_gates,"  number_active_gates(c): ",number_active_gates(c))
        println("B  number gates not equal to number active gates for goal: ",g)
        print_build_chromosome(c)
        found_c = c =  remove_inactive( c )
        num_gates = number_active_gates(found_c)
        p_current = Parameters( p.numinputs, p.numoutputs, num_gates, p.numlevelsback )
        num_gates_exceptions += 1
        complexities_list = Float64[]  # restart complexities_list
        robust_evol_list = Tuple{Float64,Float64}[]  # restart robust_evol_list 
        iter = 0
        tries = 1  # do more tries
      end
      push!(complexities_list, complexity5(c))
      push!(robust_evol_list, mutate_all( c, funcs,robustness_only=true))
      tries += 1
    end
    iter += 1
  end
  #println("tries: ",tries)
  avg_complexity = sum(complexities_list)/length(complexities_list)
  avg_robustness = sum( map(x->x[1],robust_evol_list))/length(robust_evol_list)
  avg_evolvability = sum( map(x->x[2],robust_evol_list))/length(robust_evol_list)
  (g,num_gates,number_active_gates(found_c),avg_complexity,tries,avg_robustness,avg_evolvability,num_gates_exceptions)
end
