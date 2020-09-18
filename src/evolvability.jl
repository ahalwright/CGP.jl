# Define phenotypic evolvability where a phenotype is a goal or chromosome output.
# Note that genotypic evolvability is computed by calling mutate_all( c, funcs, robust_only=true )[2]
# Use Wagner's (2008) method of evolving nchromess chromsomes whose output is the goal (where Wagner uses nchromes=100).

using DataFrames
using Dates
using CSV
using Distributed

export evolvability, run_evolvability, evo_result, test_evo, evo_result_type, run_evolve_g_pairs 
#=
mutable struct evo_result_type
  goal::Goal
  nchromes::Int64
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  maxsteps::Int64
  robustness::Int64
  evolvability::Int64
end
=#

function evo_result( g::Goal, nchromes, max_steps::Int64, p::Parameters, evolvable_count::Int64, all_count::Int64,
      evo_diff_count::Int64 )
  evo_result_type(
    g,
    nchromes,
    p.numinputs,
    p.numoutputs,
    p.numinteriors,
    p.numlevelsback,  
    max_steps,
    all_count,
    evolvable_count,
    evo_diff_count
  )
end

function evo_result( g::Goal, nchromes, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64, 
      max_steps::Int64, all_count::Int64, evolvable_count::Int64, evo_diff_count::Int64 )
  evo_result_type(
    g,
    nchromes,
    numinputs,
    numoutputs,
    numinteriors,
    numlevelsback,  
    max_steps,
    all_count,
    evolvable_count,
    evo_diff_count
  )
end

function evo_result_to_tuple( er::evo_result_type )
  (
    er.goal,
    er.nchromes,
    er.numinputs,
    er.numints,
    er.levelsback,
    er.maxsteps,
    er.all_count,
    er.evolvable_count,
    er.evo_diff_count,
    er.evolvable_count/er.all_count   # the Q field
  )
end

function evolvability( g::Goal, funcs::Vector{Func}, nchromes::Int64, maxsteps::Int64, 
      numinputs::Int64, numinteriors::Int64, numlevelsback::Int64 )
  p = Parameters( numinputs=numinputs, numoutputs=length(g), numinteriors=numinteriors, numlevelsback=numlevelsback )
  evolvability( g, funcs, nchromes,  maxsteps, p )
end

# Return dataframe of genotypic evolvability results
# If intermediate_gens is a non-empty list, save intermediate results for those number of chromosomes
function evolvability( er::evo_result_type, funcs::Vector{Func}; intermediate_gens::Vector{Int64}=Int64[] )
  repeat_limit = 10
  #println("er: ",er)
  p = Parameters( numinputs=er.numinputs, numoutputs=er.numoutputs, numinteriors=er.numints, numlevelsback=er.levelsback )
  all_sum = 0      # Accumulation of count of all neutral neighbors
  unique_result = Goal[]  # Accumulation of unique neutral neighbors
  goal_diff_counts = Int64[]  # save for intermediate generations
  if length(intermediate_gens) > 0
    @assert maximum(intermediate_gens) <= er.nchromes
    all_sums = Int64[]  # save all_sum for intermediate generations
    unique_intermeds = Vector{Goal}[]   # results saved for steps i in intermediate_gens
    intermediate_count = 1
    prev_goal_count = 0
    #println("len inter > 0")
  end
  @assert p.numoutputs==length(er.goal)
  for i = 1:er.nchromes
    numinteriors_repeat = er.numints
    numlevelsback_repeat = er.levelsback
    c = random_chromosome( p, funcs )
    (c,steps,output,mached_goals,matched_goals_list) = mut_evolve( c, [er.goal], funcs, er.maxsteps )
    repeat = 0
    while repeat < repeat_limit && steps == er.maxsteps
      numinteriors_repeat += 1
      numlevelsback_repeat += 1
      p_repeat = Parameters( numinputs=er.numinputs, numoutputs=length(er.goal), numinteriors=numinteriors_repeat, 
          numlevelsback=numlevelsback_repeat )
      println("repeating function mut_evolve with numints: ",numinteriors_repeat,"  and with levsback: ",numlevelsback_repeat)
      c = random_chromosome( p_repeat, funcs )
      (c,steps,output,mached_goals,matched_goals_list) = mut_evolve( c, [er.goal], funcs, er.maxsteps )
      repeat += 1
    end
    #print_build_chromosome( c )
    goal_list = mutate_all( c, funcs, output_outputs=true )
    all_sum += length(goal_list)
    #println("goal_list: ",goal_list)
    unique_goals = unique(goal_list)
    #print("len unique_goal_list)): ",length(unique_goals))
    unique_result = unique(vcat(unique_result,unique_goals))
    current_goal_count = length(unique_result)
    if length(intermediate_gens) > 0 && i == intermediate_gens[intermediate_count]
      diff_count = current_goal_count - prev_goal_count
      #println("  diff_count: ",diff_count)
      push!(goal_diff_counts,diff_count)
      prev_goal_count = length(unique_result)
      unique_intermed = deepcopy(unique_result)
      ##println("i: ",i,"  int_count: ",intermediate_count,"  length(intermediate): ",length(intermediate))
      push!(all_sums,all_sum)
      push!(unique_intermeds,unique_intermed)
      intermediate_count += 1
    end
  end
  if length(intermediate_gens) == 0
    #println("len(unique_result): ",length(unique_result))
    er.all_count = all_sum
    er.evolvable_count = length(unique_result)
    er.evo_diff_count = 0
    return er
  else
    er_list = evo_result_type[]
    int_count = 1
    for int in intermediate_gens
      er_i = deepcopy(er)
      er_i.nchromes = int
      er_i.all_count = all_sums[int_count]
      er_i.evolvable_count = length(unique_intermeds[int_count])
      er_i.evo_diff_count = goal_diff_counts[int_count]
      #println("int: ",int,"  er_i: ",er_i)
      push!(er_list,er_i)
      int_count += 1
    end
    return er_list
  end
end

function run_evolvability( nreps::Int64, gl::GoalList, funcs::Vector{Func}, nchromes::IntRange,  maxsteps::IntRange, 
      numinputs::Int64, numinteriors::IntRange, numlevelsback::IntRange; intermediate_gens::Vector{Int64}=Int64[] )
  if length(intermediate_gens) == 0
    evo_result_list = evo_result_type[]
  else
    evo_result_list = Vector{evo_result_type}[]
  end
  df = DataFrame()
  df.goal=Vector{MyInt}[]
  df.nchromes=Int64[]
  df.numinputs=Int64[]
  df.numints=Int64[]
  df.levsback=Int64[]
  df.maxsteps=Int64[]
  df.all_count=Int64[]
  df.evolvable_count=Int64[]
  df.evo_diff_count=Int64[]
  df.Q = Float64[]
  for g in gl
    for nch = nchromes
      for nints = numinteriors
        for levsback = numlevelsback
          for mxsteps = maxsteps
            for rep = 1:nreps
              #er = evo_result( g, nch, mxsteps, p, 0 )
              er = evo_result( g, nch, numinputs, length(g), nints, levsback, mxsteps, 0, 0, 0 )
              #er = evolvability( er, funcs )
              #push!(df,(g,nch,numinputs,nints,levsback,mxsteps,evob))
              if length(intermediate_gens)==0
                push!(evo_result_list, er )
              else
                push!(evo_result_list, [er] )
              end
            end
          end
        end
      end
    end
  end
  #println("evo_result_list: ",evo_result_list)
  if length(intermediate_gens) == 0
    new_evo_result_list = pmap( x->evolvability(x,funcs,intermediate_gens=intermediate_gens), evo_result_list )
    #new_evo_result_list = map( x->evolvability(x,funcs,intermediate_gens=intermediate_gens), evo_result_list )
  else
    new_evo_result_list = pmap( x->evolvability(x[1],funcs,intermediate_gens=intermediate_gens), evo_result_list )
    #new_evo_result_list = map( x->evolvability(x[1],funcs,intermediate_gens=intermediate_gens), evo_result_list )
  end
  #println("evo_result_list[1]: ",evo_result_list[1])
  #new_evo_result_list = [evolvability(evo_result_list[1][1],funcs,intermediate_gens=intermediate_gens)]
  for er in new_evo_result_list
    if length(intermediate_gens) == 0
      new_row = evo_result_to_tuple(er)
      push!(df, new_row)
    else
      for e in er
        new_row = evo_result_to_tuple(e)
        push!(df, new_row)
      end
    end
  end
  df
end

function run_evolvability( nreps::Int64, gl::GoalList, funcs::Vector{Func}, nchromes::IntRange,  maxsteps::IntRange, 
      numinputs::Int64, numinteriors::IntRange, numlevelsback::IntRange, csvfile::String; 
      intermediate_gens::Vector{Int64}=Int64[] )
  #result = @timed run_evolvability( nreps, gl, funcs, nchromes,  maxsteps, numinputs, numinteriors, numlevelsback ) 
  (df,ttime) = @timed run_evolvability( nreps, gl, funcs, nchromes,  maxsteps, numinputs, numinteriors, numlevelsback,
      intermediate_gens=intermediate_gens ) 
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  println("size(df): ",size(df))
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )    
  println("# date and time: ",Dates.now())
  println("# run time in minutes: ",ttime/60)
  println("# funcs: ", Main.CGP.default_funcs(numinputs[end]))
  open( csvfile, "w" ) do f 
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )    
    println(f,"# run time in minutes: ",ttime/60)
    println(f,"# funcs: ", Main.CGP.default_funcs(numinputs[end]))
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

function evolve_g_pairs( df::DataFrame, g_pair::Tuple{String,String}, p::Parameters,  maxsteps::Int64 )
  (src_g,dst_g) = g_pair
  funcs = default_funcs(p.numinputs)
  source_g = parse(MyInt,src_g)
  dest_g = parse(MyInt,dst_g)
  (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) = 
      mut_evolve_increase_numints( random_chromosome( p, funcs ), [[source_g]], funcs, maxsteps )
  if step == maxsteps
    println("failed to find chromosome for source goal: ",src_g)
    return nothing
  else
    println("found source goal: ",src_g)
  end
  (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) = 
      mut_evolve_increase_numints(c, [[dest_g]], funcs, maxsteps )
  if step == maxsteps
    println("failed to find chromosome for dest goal: ",dst_g)
    return nothing
  else
    println("found dest goal: ",dst_g)
  end
  total_steps = (new_levsback - p.numlevelsback)*maxsteps + step
  s_df = df[df.goal.==src_g,[:counts10ints,:complex,:mutrobust,:evolvability]]
  d_df = df[df.goal.==dst_g,[:counts10ints,:complex,:mutrobust,:evolvability]]
  hdist = hamming_distance(source_g,dest_g,p.numinputs)
  # changed order on 9/15/20
  row=(src_g,s_df[1,:counts10ints],s_df[1,:complex],s_df[1,:mutrobust],s_df[1,:evolvability],
        dst_g,d_df[1,:counts10ints],d_df[1,:complex],d_df[1,:mutrobust],d_df[1,:evolvability],hdist,
        p.numinputs,p.numoutputs,new_numints,new_levsback,maxsteps,total_steps)
  row
end

# Test Wagner's hypothesis that "An evolutionary search’s ability to find a target genotype is only weakly 
#    correlated with the evolvability of the source genotype".
# Choose a random sample of sample_genotypes for source and destination genotypes.
# For nreps pairs (source_g, dest_g) from this sample, calulate the genotypic evolvability if not already done.
# Then for each pair run mut_evolve from a genotype that corresponds to the source_g with goal dest_g nruns times,
#   and count the number of steps for each run.  Use a large maxsteps and a large numints so that
#   reruns due to failures are rare.
# Returns a DataFrame whose fields are given below.  
# Objective: determine correlation of steps with src_count, src_complex, dst_count, dst_complex
# Wagner (2008) claims that there is little correlation of difficulty with src for his RNA data
# As of 9/15/20, the dataframe df can be obtained by 
# df = read_dataframe("../data/consolidate/geno_pheno_raman_df_all_9_13.csv");
function run_evolve_g_pairs( df::DataFrame, sample_size::Int64, nreps::Int64, nruns::Int64, numints::Int64, maxsteps::Int64 )
  println("run_evolve_g_pairs")
  p = Parameters( numinputs=df.numinputs[1], numoutputs=df.numoutputs[1], numinteriors=numints, numlevelsback=df.levsback[1] )
  funcs = default_funcs(p.numinputs)
  ndf = DataFrame()
  ndf.source_g = String[]
  ndf.src_count=Int64[]
  ndf.src_cmplx=Float64[]
  ndf.src_mrobust=Float64[]
  ndf.src_evolble=Float64[]
  ndf.dest_g = String[]
  ndf.dst_count=Int64[]
  ndf.dst_cmplx=Float64[]
  ndf.dst_mrobust=Float64[]
  ndf.dst_evolble=Float64[]
  ndf.hamming_dist = Float64[]
  ndf.numinputs = Int64[]
  ndf.numoutputs = Int64[]
  ndf.numints = Int64[]
  ndf.levelsback = Int64[]
  ndf.maxsteps = Int64[]
  ndf.steps = Int64[]
  samples = rand( 1:size(df)[1], sample_size )
  row_list = Tuple[]
  g_pair_list = Tuple{String,String}[]
  for i = 1:nreps
    src_g = df.goal[rand(samples)]
    dst_g = df.goal[rand(samples)]
    for j = 1:nruns
      push!(g_pair_list,(src_g,dst_g))
    end
  end
  row_list = pmap( pair->evolve_g_pairs( df, pair, p, maxsteps), g_pair_list )
  #row_list = map( pair->evolve_g_pairs( df, pair, p, maxsteps), g_pair_list )
  for row in row_list
    if row != nothing
      push!(ndf,row)
    else
      println("evolve_g_pairs() failed!")
    end
  end
  ndf
end

function run_sample_g_pairs( df::DataFrame, nreps::Int64, sample_size::Int64, order_by::Symbol,
    p::Parameters, maxsteps::Int64 )
  pmap( x-> sample_g_pairs( df, nreps, sample_size, order_by, p, maxsteps, x ), [1,2,3,4] )
end 

# Choose a random sample of genotypes of size 4*sample size.  
# Sort this sample according to field order_by of DataFrame df.  Then return the smallest eighth
#  and the largest eigth
function sample_g_pairs( df::DataFrame, nreps::Int64, sample_size::Int64, order_by::Symbol, 
    p::Parameters, maxsteps::Int64, option::Int64 )
  sample_size_multiplier = 8
  sdf = DataFrame()
  sdf.goals = rand(MyInt(0):MyInt(size(df)[1]-1),sample_size_multiplier*sample_size)
  #sdf.goals = map(x->@sprintf("0x%x",x), rand(1:size(df)[1],sample_size_multiplier*sample_size))
  #println("sdf.goals: ",sdf.goals)
  sdf[!,order_by] = [ df[ gg, order_by ] for gg in randgoallist(sample_size_multiplier*sample_size,p.numinputs,p.numoutputs) ] 
  nsdf = sort(sdf,order(order_by))
  simple_goals = nsdf.goals[1:sample_size]
  complex_goals = nsdf.goals[((sample_size_multiplier-1)*sample_size+1):(sample_size_multiplier*sample_size)]
  steps_list = Int64[]
  steps_count = 0
  steps_squared_count = 0
  for i = 1:nreps
    if option == 1
      src_g = rand(simple_goals)
      dst_g = rand(simple_goals)
    elseif option == 2
      src_g = rand(simple_goals)
      dst_g = rand(complex_goals)
    elseif option == 3
      src_g = rand(complex_goals)
      dst_g = rand(simple_goals)
    elseif option == 4
      src_g = rand(simple_goals)
      dst_g = rand(complex_goals)
    end
    steps = evolve_g_pair( [src_g], [dst_g], p, maxsteps)
    #println("steps: ",steps,"  steps_count: ",steps_count)
    push!(steps_list,steps)
    steps_count += steps
    steps_squared_count += steps^2
  end
  variance = (steps_squared_count - steps_count^2/nreps)/(nreps-1)
  (option, steps_count/nreps, variance, steps_list )
end

function evolve_g_pair( s_g::Goal, d_g::Goal, p::Parameters, maxsteps::Int64 )
  funcs = default_funcs(p.numinputs)
  c = random_chromosome( p, funcs)
  (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) =
      mut_evolve_increase_numints(c, [s_g], funcs, maxsteps ) 
  println("src crhomsome evolved in ",step," steps.")
  (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) =
      mut_evolve_increase_numints(c, [d_g], funcs, maxsteps ) 
  println("dst crhomsome evolved in ",step," steps.")
  step
end

function run_evolve_g_pairs( df::DataFrame, sample_size::Int64, nreps::Int64, nruns::Int64, 
      numints::Int64, maxsteps::Int64, csvfile::String)
  println("run_evolve_g_pairs: csvfile: ",csvfile)
  p = Parameters( numinputs=df.numinputs[1], numoutputs=df.numoutputs[1], numinteriors=numints, numlevelsback=df.levsback[1] )
  (ndf,ttime) = @timed run_evolve_g_pairs( df, sample_size, nreps, nruns, numints, maxsteps )  
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  println("size(df): ",size(df))
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )    
  println("# date and time: ",Dates.now())
  println("# run time in minutes: ",ttime/60)
  println("# funcs: ", Main.CGP.default_funcs(p.numinputs))
  open( csvfile, "w" ) do f 
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )    
    println(f,"# run time in minutes: ",ttime/60)
    println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
    CSV.write( f, ndf, append=true, writeheader=true )
  end
  ndf
end

# Test the Wagner (2008) that genotypic robustness is negatively associated with genotypic evolvability.  
# Strongly confirmed by:
# julia> tge=test_g_evolvability_g_robustness( 100000, 4, 1, 9, 5 );
# julia> corspearman( tge[1], tge[2])
#  -0.8104812361089767
function test_g_evolvability_g_robustness( nreps::Int64, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64 )
  p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
  funcs = default_funcs(numinputs)
  robust_vec = Float64[]
  evolvable_vec = Float64[]
  for i = 1:nreps
    c = random_chromosome(p,funcs)
    c_output = output_values(c)
    #println("c_output: ",c_output)
    outputs = mutate_all( c, funcs, output_outputs=true )
    #println("outputs: ",outputs)
    robust_outputs = filter( x->x==c_output, outputs )
    #println("robust_outputs: ",robust_outputs)
    evolvable_outputs = unique(outputs)
    #println("evolvable_outputs: ",evolvable_outputs)
    #println("push: ",length(robust_outputs)/length(outputs),"  ",length(evolvable_outputs)/length(outputs))
    push!( robust_vec, length(robust_outputs)/length(outputs))
    push!( evolvable_vec, length(evolvable_outputs)/length(outputs))
  end
  ( robust_vec, evolvable_vec )
end

function run_geno_robustness( ngoals::Int64, nreps::Int64, numinputs::Int64, numoutputs::Int64, 
    numinteriors::Int64, numlevelsback::Int64, max_steps::Int64 )
  robust_vec = Float64[]
  evolvable_vec = Float64[]
  for j = 1:ngoals
    goal = randgoal( numinputs, numoutputs)
    println("goal: ",goal)
    (rbst, evbl) = geno_robustness( goal, nreps, numinputs, numoutputs, numinteriors, 
        numlevelsback, max_steps )
    push!(robust_vec,rbst)
    push!(evolvable_vec,evbl)
  end
  (robust_vec, evolvable_vec)
end

function geno_robustness( goal::Goal, nreps::Int64, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64, max_steps::Int64 )
  p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
  funcs = default_funcs(numinputs)
  robust_vec = Float64[]
  evolvable_vec = Float64[]
  for i = 1:nreps
    c = random_chromosome(p,funcs)
    c_output = output_values(c)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) =
      mut_evolve_increase_numints( c, [goal], funcs, max_steps )
    c_output = output_values(c)
    @assert c_output == goal
    #println("output values: ",c_output)
    outputs = mutate_all( c, funcs, output_outputs=true )
    robust_outputs = filter( x->x==c_output, outputs )
    evolvable_outputs = unique(outputs)
    push!( robust_vec, length(robust_outputs)/length(outputs))
    push!( evolvable_vec, length(evolvable_outputs)/length(outputs))
  end
  return ( sum(robust_vec)/length(robust_vec), sum(evolvable_vec)/length(evolvable_vec) )
end

function test_evo()
  g = [0x0332]
  numinputs=4; numoutputs=1; numinteriors=9; numlevelsback=5;
  setup_funcs(numinputs)
  g =[0x03f3]
  p = Parameters( numinputs=numinputs, numoutputs=length(g), numinteriors=numinteriors, numlevelsback=numlevelsback )
  max_steps = 200000
  nchromes = 50
  evo_result = Main.CGP.evo_result
  er = evo_result( g, nchromes, max_steps, p, 0, 0, 0 )
  funcs=default_funcs(p.numinputs)
  #evolvability( er, funcs )
  #evolvability( er, funcs, intermediate_gens=[1,4,8] )
  run_evolvability( 2, [g], funcs, er.nchromes, er.maxsteps, er.numinputs, er.numints, er.levelsback, intermediate_gens=[10,20,30,40,50] )
  #run_evolvability( 20, [g], funcs, 5:5:20, er.maxsteps, er.numinputs, er.numints, er.levelsback, "test.csv" )
end

