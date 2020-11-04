# Includes explore_complexity()
using DataFrames, Statistics, Dates, CSV, Distributed

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
        quantile(complexity_list,0.90)
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
      csvfile::String="" )
  goallist = randgoallist( ngoals, p.numinputs,p.numoutputs)
  run_explore_complexity( nreps, nruns, p, goallist, num_circuits, max_ev_steps, csvfile=csvfile )
end

# Run the next version of run_explore_complexity() nreps time parallelized by pmap().
# This version takes goallist as an argument.
function run_explore_complexity( nreps::Int64, nruns::Int64, p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64;
      csvfile::String="" )
  result_dfs = DataFrame[]
  result_dfs = pmap(x->run_explore_complexity( nruns, p, goallist, num_circuits, max_ev_steps), collect(1:nreps))
  #result_dfs = map(x->run_explore_complexity( nruns, p, goallist, num_circuits, max_ev_steps), collect(1:nreps))
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
function run_explore_complexity( nruns::Int64, p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64)
  df = DataFrame()
  df.steps = Int64[]
  df.unique_goals = Int64[]
  df.mean_complexity = Float64[]
  df.max_complexity = Float64[]
  (circuit_list,goals_found,row_tuple) = explore_complexity( p, goallist, num_circuits, max_ev_steps )
  push!(df,row_tuple)
  reduced_goallist = goallist
  for i = 2:(nruns-1)
    extended_circuit_list = extend_list_by_dups( circuit_list, num_circuits )
    reduced_goallist = setdiff(reduced_goallist,goals_found)
    #println("run i: ",i,"  length(circuit_list): ",length(circuit_list),"  length(reduced_goallist): ",length(reduced_goallist))
    (circuit_list,goals_found,row_tuple) = explore_complexity( p, reduced_goallist, extended_circuit_list, max_ev_steps )
    push!(df,row_tuple)
  end
  extended_circuit_list = extend_list_by_dups( circuit_list, num_circuits )
  reduced_goallist = setdiff(reduced_goallist,goals_found)
  (circuit_list,goals_found,row_tuple) = explore_complexity( p, reduced_goallist, extended_circuit_list, max_ev_steps )
  println("explore with nruns: ",nruns) 
  println("row_tuple: ",row_tuple)
  push!(df,row_tuple)
  df
end

function explore_complexity( p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64 )
  funcs = default_funcs(p.numinputs)
  circuit_list = [ random_chromosome(p,funcs) for _=1:num_circuits ]
  explore_complexity( p, goallist, circuit_list, max_ev_steps )
end

function explore_complexity( p::Parameters, goallist::GoalList, num_circuits::Int64, max_ev_steps::Int64 )
  funcs = default_funcs(p.numinputs)
  circuit_list = [ random_chromosome(p,funcs) for _=1:num_circuits ]
  explore_complexity( p, goallist, circuit_list, max_ev_steps )
end
  
# This is the explore_complexity() function that runs one generation starting with ciruits in circuit_list.
# max_ev_steps  is the maximum number of steps in the evolution of a circuit to any of the goals in goalllist.
# Returns a triple (circuits_list, goals_found, row_tuple) where:
#   circuits_list is the list of evolved circuits that output one of the goals
#   goals_found is the list of goals that are found by evolution
#   row_tuple is the row of statistics that can be pushed onto a dataframe
function explore_complexity( p::Parameters, goallist::GoalList, circuit_list::Vector{Chromosome}, max_ev_steps::Int64 )
  funcs = default_funcs(p.numinputs)
  goals_found = Goal[]
  steps_list = Int64[]
  src_cmplx_list = Float64[]
  dest_cmplx_list = Float64[]
  circuits_list = Chromosome[]
  #if length(circuit_list) > 0
    for c in circuit_list
      (new_c,step,worse,same,better,output,matched_goals,matched_goals_list) = 
          mut_evolve( c, goallist, funcs, max_ev_steps, print_steps=false )
      if step < max_ev_steps
        #print("i: ",i,"  steps: ",step)
        #@printf("  goal: 0x%04x\n",matched_goals_list[1][1][1])
        push!(goals_found, [matched_goals_list[1][1][1]])
        push!(steps_list,step)
        push!(src_cmplx_list, complexity5(c))
        push!(dest_cmplx_list, complexity5(new_c))
        push!(circuits_list, new_c )
      else
        push!(goals_found, Goal[] )
        push!(steps_list,step)
        push!(src_cmplx_list, complexity5(c))
        push!(dest_cmplx_list, 0.0)
        #push!(circuits_list, Chromosome[] )
      end
    end
  #end
  if length(steps_list) > 0
    row_tuple = (length(steps_list), length(unique(goals_found)),mean(dest_cmplx_list), maximum(dest_cmplx_list) )
  else
    row_tuple = (0,0,0.0,0.0)
  end
  println("length(circuits_list: ",length(circuits_list))
  return (circuits_list, goals_found, row_tuple)
end

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
  
