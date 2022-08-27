using DataFrames, Printf, Distributed, Dates, CSV
# Choose a random sample of goals and evolve to find one of them.
# Use two sets of goals from dataframe src_df and dest_df (whose properties have already been computed).
# Evolve from each goal a random sample goal of size sample_size from src_df
# to each goal in dest_df.
# Compare the properties of the source and destination goals.
# n_repeats is number of tries of mut_evolve() when step limit is reached.
# Suitable dataframe files:  "../data/10_27/geno_complexity10_27Fncons.csv" and "../data/10_27/geno_complexity10_27Gccons.csv"
function run_evolve_to_random_goals( p::Parameters, dest_df::DataFrame, src_df::DataFrame, sample_size::Int64, maxsteps::Int64,
    dest_csv::String="", src_csv::String=""; csvfile::String="", n_repeats::Int64=20 )
  #n_repeats = 20  # number of tries of mut_evolve().  Could make this a keyword parameter.
  funcs = default_funcs(p.numinputs)
  if size(dest_df)[1] == 0
    dest_df = read_dataframe( dest_csv )
  end
  if size(src_df)[1] == 0
    src_df = read_dataframe( src_csv )
  end  
  src_indices = rand(1:size(src_df)[1],sample_size)
  src_goallist = map(x->eval(Meta.parse(x)),src_df.goal[src_indices])   # Convert from string format to julia expression
  if src_df == dest_df
    dest_indices = setdiff(collect(1:size(src_df)[1]),src_indices)  # Do not allow duplicate source and dest goals
  else
    dest_indices = collect(1:size(srd_df)[1])
  end
  println("length(goal_indices): ",length(dest_indices))
  dest_goallist = map(x->eval(Meta.parse(x)),dest_df.goal[dest_indices])   # Convert from string format to julia expression
  println("intersect(src_goallist,dest_goallist): ",intersect(src_goallist,dest_goallist))
  @assert length(intersect(src_goallist,dest_goallist)) == 0
  df = DataFrame()
  df.src_goal = Goal[]
  df.dest_goal = Goal[]
  df.steps = Int64[]
  df.hamming = Float64[]
  df.src_complexity = Float64[]
  df.dest_complexity = Float64[]
  df.src_evo_count = Float64[]
  df.dest_evo_count = Float64[]
  df.src_robustness = Float64[]
  df.dest_robustness = Float64[]
  df.src_freq= Int64[]
  df.log_src_freq= Float64[]
  df.dest_freq = Int64[]
  df.log_dest_freq= Float64[]
  row_tuple_list = pmap( x->evolve_to_random_goal( p, dest_goallist, src_goallist, dest_df, src_df, sample_size, maxsteps, n_repeats=n_repeats), collect(1:sample_size))
  #row_tuple_list = map( x->evolve_to_random_goal( p, dest_goallist, src_goallist, dest_df, src_df, sample_size, maxsteps, n_repeats=n_repeats), collect(1:sample_size))
  for row_tuple in row_tuple_list
    #println("size(df): ",size(df),"  length(row_tuple): ",length(row_tuple))
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
      CSV.write( f, df, append=true, writeheader=true ) 
    end
  else
    println("# date and time: ",Dates.now())
    println("# host: ",hostname," with ",nprocs()-1,"  processes: " )   
    println("# funcs: ", Main.CGP.default_funcs(p.numinputs))
    print_parameters(p,comment=true)
    println("# maxsteps: ",maxsteps)      
    df
  end
end

function evolve_to_random_goal( p::Parameters, dest_goallist::GoalList, src_goallist::GoalList, dest_df::DataFrame, src_df::DataFrame, sample_size::Int64, maxsteps::Int64;
      n_repeats::Int64=20 )
  funcs = default_funcs(p.numinputs)
  c = random_chromosome(p,funcs)
  # Find a circuit that maps to g
  i = rand(1:length(src_goallist))
  #println("[src_goallist[i]]: ",[src_goallist[i]])
  (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve_repeat( n_repeats, p, [src_goallist[i]], funcs, maxsteps )
  #print_build_chromosome(c)
  # Find a circuit that maps to one of the goals in dest_goallist
  #println("dest_goallist[1:5]: ",dest_goallist[1:5])
  #(c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve_repeat( n_repeats, p, dest_goallist, funcs, maxsteps )
  count = 0
  total_steps = 0
  (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, dest_goallist, funcs, maxsteps )
  total_steps += step
  while step == maxsteps && count < n_repeats
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, dest_goallist, funcs, maxsteps )
    total_steps += step
    count += 1
  end
  if step == maxsteps || count == n_repeats
    println("step == maxsteps in function evolve_to_random_goals()")  # this message doesn't print, instead get message: no method matching iterate(::Nothing)
  end
  dest_goal = [matched_goals_list[1][1][1]]
  #println("type dest_goal: ",typeof(dest_goal),"  dest_goal: ",dest_goal)
  dest_goal_str = @sprintf("UInt16[0x%04x]",dest_goal[1])
  #println("dest_g: ",dest_goal_str)
  hamming_dist = hamming_distance( src_goallist[i], dest_goal, p.numinputs )
  src_complexity = src_df[i,:complexity]
  src_evo_count = src_df[i,:evo_count]
  src_robustness = src_df[i,:robustness]
  #println("src: ni: ",src_df.numints[i],"  lb: ",src_df.numlevsback[i],"  nm11: ",:ints11_8 in names(src_df),"  nm10: ",:ints10_5 in names(src_df))
  if src_df.numints[i]==11 && src_df.numlevsback[i]==8 && :ints11_8 in names(src_df)
    src_freq = src_df.ints11_8[i]
  elseif src_df.numints[i]==10 && src_df.numlevsback[i]==5 && :ints10_5 in names(src_df)
    src_freq = src_df.ints10_5[i]
  else
    src_freq = 0.0
  end
  log_src_freq = src_freq > 0 ? log10(src_freq) : 0.0
  if dest_df.numints[1]==11 && dest_df.numlevsback[1]==8 && :ints11_8 in names(dest_df)
    dest_freq = dest_df[dest_df.goal.==dest_goal_str,:ints11_8][1]
  elseif dest_df.numints[1]==10 && dest_df.numlevsback[1]==5 && :ints10_5 in names(dest_df)
    dest_freq = dest_df[dest_df.goal.==dest_goal_str,:ints10_5][1]
  else
    dest_freq = 0.0
  end
  log_dest_freq = dest_freq > 0 ? log10(dest_freq) : 0.0
  dest_complexity = dest_df[dest_df.goal.==@sprintf("UInt16[0x%04x]",dest_goal[1]),:complexity][1]
  dest_evo_count = dest_df[dest_df.goal.==@sprintf("UInt16[0x%04x]",dest_goal[1]),:evo_count][1]
  dest_robustness = dest_df[dest_df.goal.==@sprintf("UInt16[0x%04x]",dest_goal[1]),:robustness][1]
  row_tuple = 
  ( 
    src_goallist[i], 
    dest_goal, 
    total_steps, 
    hamming_dist,
    src_complexity, 
    dest_complexity,
    src_evo_count, 
    dest_evo_count,
    src_robustness, 
    dest_robustness,
    src_freq,
    log_src_freq,
    dest_freq,
    log_dest_freq
  )
  #println("row_tuple: ",row_tuple)
  row_tuple
end

#= Example call:
p = Parameters(4,1,11,8) 
sdf = read_dataframe("../data/10_27/geno_complexity10_27FMNccons.csv")
ddf = read_dataframe("../data/10_27/geno_complexity10_27FMNccons.csv")
run_evolve_to_random_goals( p, ddf, sdf, 8, 100000)

300434033794290
06965
=#
