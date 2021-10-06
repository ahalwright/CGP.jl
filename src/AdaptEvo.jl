# Project on adaptive evolvability as first described in cce/Adaptive_evolvability9_6_21.docx.

# For each of number_sources starting phenotypes choose number_target random phenotypes within 
#     distance_limit of the source phenotype.
# For each starting phenotype and each target phenotype do number_runs epochal evolutions with 
#     a maximum of step_limit steps and fitness defined by hamming distance from target.
# Record the number of successes and the averge number of steps for successful runs
# max_steps is the maximum number of steps on the evolution to find the starting chromosome for a run
function compare_evolvabilities( p::Parameters, source_list::GoalList, number_targets::Int64, 
    distance_limit::Float64, step_limit::Int64, number_runs::Int64;
    csvfile::String="", evo_df_file::String="", evo_df_column::Symbol=:d_evolvability, 
    iteration_limit::Int64=200, max_steps::Int64=200000, distance_function::Function=hamming_distance )
  evo_df = DataFrame()
  if length(evo_df_file) > 0
    try
      evo_df = read_dataframe(evo_df_file)
    catch
      error("reading dataframe from ",evo_df_file," failed")
    end
    @assert size(evo_df)[1] == length(source_list)  # Make sure goal lists are compatible
  end
  df = DataFrame()
  df.source = Goal[]
  #df.numinputs = Int64[]
  #df.numoutputs = Int64[]
  #df.numgates= Int64[]
  #df.numlevelsback = Int64[]
  df.numtargets = Int64[]
  df.dist_limit = Float64[]
  df.step_limit = Int64[]
  df.number_runs = Int64[]
  df.avg_num_successes = Float64[]
  df.avg_steps = Float64[]
  println("size(df): ",size(df))
  # list of parameter values that go into the dataframe
  #df_params_list = Any[p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, number_targets, distance_limit, step_limit, number_runs]
  #println("df_params_list: ",df_params_list)
  #for src in source_list
  df_rows = pmap( src->run_compare_evolvabilities( p, src, number_targets, distance_limit, step_limit, number_runs,
        iteration_limit=iteration_limit, max_steps=max_steps, distance_function=distance_function ), source_list )
  for df_row in df_rows
    push!(df,df_row)
  end
  df.evolvability = evo_df[:,evo_df_column]
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Helper function for the above function for a single src goal.
function run_compare_evolvabilities( p::Parameters, src::Goal, number_targets::Int64, 
  distance_limit::Float64, step_limit::Int64, number_runs::Int64;
  iteration_limit::Int64=200, max_steps::Int64=200000, distance_function::Function=hamming_distance )
  println("src: ",src)
  df_params_list = 
    Any[
      #p.numinputs,
      #p.numoutputs, 
      #p.numinteriors, 
      #p.numlevelsback, 
      number_targets, 
      distance_limit, 
      step_limit, 
      number_runs
    ]
  iteration = 0
  target_list = Goal[]
  while iteration < iteration_limit && length(target_list) < number_targets
    target = randgoal( p.numinputs, p.numoutputs )
    #hd = hamming_distance( src, target, p.numinputs ) 
    hd = distance_function( src, target, p.numinputs ) 
    #println("target: ",target,"  hd: ",hd)
    if hd <= distance_limit && target != src
      push!(target_list,target)
      println("target: ",target,"  hd: ",hd,"  iteration: ",iteration)
    end
    iteration += 1
  end
  if iteration == iteration_limit
    error("target list iteration exceeded iteration limit of ",iteration_limit," for source src: ",src)
  end
  println("target_list: ",target_list)
  successes = 0
  sum_steps = 0
  for target in target_list
    ch = random_chromosome( p )
    (src_ch,steps) = neutral_evolution(ch, src, max_steps )
    while iteration < iteration_limit && steps == max_steps
      (src_ch,steps) = neutral_evolution(ch, src, max_steps )
      iteration += 1
    end
    if iteration == iteration_limit
      error("iteration to src exceeded iteration limit of ",iteration_limit," for source src: ",src," and target: ",target)
    end
    @assert  steps < max_steps
    for run = 1:number_runs
      (target_ch,steps) = neutral_evolution(src_ch, target, step_limit )        
      successes += steps < step_limit ? 1 : 0
      sum_steps += steps < step_limit ? steps : 0
    end
  end
  avg_successes = successes/length(target_list)/number_runs
  avg_steps = sum_steps/successes
  df_row = vcat([src],df_params_list,[avg_successes,avg_steps])
  return df_row
end

# Doesn't make these assignments global
function test_inputs()
  p = Parameters(3,1,7,4)
   #source_list=randgoallist( 10, p.numinputs, p.numoutputs )
   source_list=map(x->[MyInt(x)],collect(0:255))
   number_targets=3
   distance_limit=0.3
   step_limit=200
   number_runs=10
   evo_df_file="../data/9_3_21/robust_evo_by_walks9_3_21B.csv"
end

