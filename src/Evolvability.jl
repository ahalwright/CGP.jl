# Define phenotypic evolvability where a phenotype is a goal or chromosome output.
# Note that genotypic evolvability is computed by calling mutate_all( c, funcs, robust_only=true )[2]
# Use Wagner's (2008) method of evolving nchromess chromsomes whose output is the goal (where Wagner uses nchromes=100).

using DataFrames
using Dates
using CSV
using Distributed
using HypothesisTests
using Statistics
using Printf

export evolvability, run_evolvability, evo_result, test_evo, evo_result_type, run_evolve_g_pairs 
export run_geno_robustness, geno_robustness, evo_robust, random_neutral_walk, run_random_neutral_walk
export run_geno_complexity, geno_complexity, pop_evolvability_robustness
export parent_child_complexity, rand_norm, scatter_plot
export pheno_set_rand_neutral_walks, pheno_set_rand_neutral_walk
#=  Moved to aliases.jl so that this file can be included.
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

function evo_result( g::Goal, nchromes, max_steps::Int64, p::Parameters, nrepeats::Int64, evolvable_count::Int64, all_count::Int64,
      evo_diff_count::Int64 )
  evo_result_type(
    g,
    nchromes,
    p.numinputs,
    p.numoutputs,
    p.numinteriors,
    p.numlevelsback,  
    nrepeats,
    max_steps,
    all_count,
    evolvable_count,
    evo_diff_count
  )
end

function evo_result( g::Goal, nchromes, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64, 
      max_steps::Int64, nrepeats::Int64, all_count::Int64, evolvable_count::Int64, evo_diff_count::Int64 )
  evo_result_type(
    g,
    nchromes,
    numinputs,
    numoutputs,
    numinteriors,
    numlevelsback,  
    nrepeats,
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
    er.nrepeats,
    er.all_count,
    er.evolvable_count,
    er.evo_diff_count,
    er.evolvable_count/er.all_count   # the Q field
  )
end

# return the pair (pheno_evolvability, pheno_robustness)
function evo_robust( g::Goal, p::Parameters, funcs::Vector{Func}, nchromes::Int64, maxsteps::Int64, n_repeats::Int64  )
  all_count = 0
  robust_count = 0
  goal_list = Goal[]
  for i = 1:nchromes
    j = 1
    while j < n_repeats
      res =  mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps)
      if res == nothing
        j += 1
        continue
      end
      (c,step,worse,same,better,output,matched_goals,matched_goals_list) = res
      @assert output_values(c) == g
      #print("i: ",i,"  found c: ")
      #print_build_chromosome(c)
      outputs = mutate_all( c, funcs, output_outputs=true )
      #println("outputs: ",outputs)
      #println("i: ",i,"  unique outputs: ",sort(unique(outputs)))
      all_count += length(outputs)
      robust_count += length(filter( x->x == g, outputs ))
      #println("len unique outputs: ", length(unique(outputs)))
      #println("robust outputs: ",length(filter( x->x == g, outputs )))
      goal_list = unique(vcat( goal_list, unique(outputs) ))
      #println("goal_list: ",sort(goal_list))
      println((all_count, robust_count, length(goal_list )))
      j = n_repeats
    end
  end
  println((all_count, robust_count, length(goal_list )))
  ( length(goal_list)/all_count, robust_count/all_count )
end
    
function evolvability( g::Goal, funcs::Vector{Func}, nchromes::Int64, maxsteps::Int64, 
      numinputs::Int64, numinteriors::Int64, numlevelsback::Int64 )
  p = Parameters( numinputs=numinputs, numoutputs=length(g), numinteriors=numinteriors, numlevelsback=numlevelsback )
  evolvability( g, funcs, nchromes,  maxsteps, p )
end

# Update er with results of one run of computation of evolvability
# The following is done nchromes times:
#   Start with a random chromosome and evolve the goal.
#   When the goal is found, generate the output of the chromosome.
# Evolvability is the count of the unique goals over all of the nchromes chromosome outputs.
# TODO:  convert this to a relative evolvability.  
# If intermediate_gens is a non-empty list, save intermediate results for those number of chromosomes
function evolvability( er::evo_result_type, funcs::Vector{Func}, max_repeats::Int64; intermediate_gens::Vector{Int64}=Int64[],
    increase_numints::Bool=true )
  #println("er: ",er)
  succeed = true    # Whether the evoluton succeeded.  Reset below
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
  er.nrepeats = 0    # Establish scope for variable nrepeats
  @assert p.numoutputs==length(er.goal)
  for i = 1:er.nchromes
    if increase_numints
      numinteriors_repeat = er.numints
      numlevelsback_repeat = er.levelsback
    end
    c = random_chromosome( p, funcs )
    (c,steps,output,mached_goals,matched_goals_list) = mut_evolve( c, [er.goal], funcs, er.maxsteps )
    nrepeats = 0
    while nrepeats < max_repeats && steps == er.maxsteps
      if increase_numints
        numinteriors_repeat += 1
        numlevelsback_repeat += 1
        p_repeat = Parameters( numinputs=er.numinputs, numoutputs=length(er.goal), numinteriors=numinteriors_repeat, 
            numlevelsback=numlevelsback_repeat )
        println("repeating function mut_evolve with numints: ",numinteriors_repeat,"  and with levsback: ",numlevelsback_repeat)
      else
        p_repeat = p
        println("repeating function mut_evolve" )
      end
      c = random_chromosome( p_repeat, funcs )
      (c,steps,output,mached_goals,matched_goals_list) = mut_evolve( c, [er.goal], funcs, er.maxsteps )
      nrepeats += 1
    end
    succeed = (steps < er.maxsteps)
    if !succeed
      println("nrepeats: ",nrepeats)
      println("evolution failed for goal: ",er.goal)
      er.nrepeats += 100000
    end
    #print_build_chromosome( c )
    goal_list = succeed ? mutate_all( c, funcs, output_outputs=true ) : Vector{MyInt}[] # goals for all possible mutations
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
    er.nrepeats += nrepeats
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
      numinputs::Int64, numinteriors::IntRange, numlevelsback::IntRange, max_repeats::Int64; 
      increase_numints::Bool=false, intermediate_gens::Vector{Int64}=Int64[] )
  #repeats=3
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
  df.nrepeats=Int64[]
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
              er = evo_result( g, nch, numinputs, length(g), nints, levsback, mxsteps, 0, 0, 0, 0 )
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
    #new_evo_result_list = pmap( x->evolvability(x,funcs,max_repeats,increase_numints=increase_numints,intermediate_gens=intermediate_gens), evo_result_list )
    new_evo_result_list = map( x->evolvability(x,funcs,max_repeats,increase_numints=increase_numints,intermediate_gens=intermediate_gens), evo_result_list )
  else
    #new_evo_result_list = pmap( x->evolvability(x[1],funcs,max_repeats,increase_numints=increase_numints,intermediate_gens=intermediate_gens), evo_result_list )
    new_evo_result_list = map( x->evolvability(x[1],funcs,max_repeats,increase_numints=increase_numints,intermediate_gens=intermediate_gens), evo_result_list )
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

function run_evolvability( nreps::Int64, gl::GoalList, funcs::Vector{Func}, nchromes::IntRange, maxsteps::IntRange, 
      numinputs::Int64, numinteriors::IntRange, numlevelsback::IntRange, max_repeats::Int64, csvfile::String; 
      increase_numints::Bool=false, intermediate_gens::Vector{Int64}=Int64[] )
  #result = @timed run_evolvability( nreps, gl, funcs, nchromes,  maxsteps, numinputs, numinteriors, numlevelsback ) 
  (df,ttime) = @timed run_evolvability( nreps, gl, funcs, nchromes,  maxsteps, numinputs, numinteriors, numlevelsback,
      max_repeats, intermediate_gens=intermediate_gens ) 
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  println("size(df): ",size(df))
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )    
  println("# date and time: ",Dates.now())
  println("# run time in minutes: ",ttime/60)
  println("# max_repeats: ",max_repeats)
  println("# funcs: ", Main.CGP.default_funcs(numinputs[end]))
  open( csvfile, "w" ) do f 
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )    
    println(f,"# run time in minutes: ",ttime/60)
    println(f,"# max_repeats: ",max_repeats)
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

# Calls the above version of run_evolve_g_pairs() with a csvfile argument.
# Outputs the resulting dataframe to the csvfile.
# For example runs, see data/9_15, 9_18, 10_8, 10_29
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

# Test hypotheses that evolution starting with complex goals takes fewer cases than evolution starting with simple goals.
# There are two cases:  evolution of simple goals and evolution of complex goals.
# df is the dataframe that results from running 10 (in the 4x1) case evolutions to find each goal.
# df = read_dataframe("../data/consolidate/geno_pheno_raman_df_all_9_13.csv");
# df = read_dataframe("../../../complexity/data/consolidate/geno_pheno_raman_df_all_9_13.csv");
# sample_size is the size of the sample used in this run.  See the comments on sample_g_pairs()
# order_by is the column of the dataframe df used to measure complexity.  Default is :complex
# Defaults for the 4x1 case might be numinteriors=10, numlevelsback=5, maxsteps=150000
function run_sample_g_pairs( df::DataFrame, nreps::Int64, sample_size::Int64, order_by::Symbol,
    p::Parameters, maxsteps::Int64; hamming::Bool=false )
  result = pmap( x-> sample_g_pairs( df, nreps, sample_size, order_by, p, maxsteps, x, hamming=hamming ), [1,2,3,4] )
  #result = map( x-> sample_g_pairs( df, nreps, sample_size, order_by, p, maxsteps, x, hamming=hamming ), [1,2,3,4] )
  pvalues = 
      (pvalue(EqualVarianceTTest(result[1][4],result[3][4])),
       pvalue(EqualVarianceTTest(result[2][4],result[4][4])))
  (pvalues,result)
end 

# Choose a random sample of genotypes of size 4*sample size.  
# Sort this sample according to field order_by of DataFrame df.  Then return the smallest eighth
#  and the largest eigth (assuming that sample_size_multiplier==8 as set below).
function sample_g_pairs( df::DataFrame, nreps::Int64, sample_size::Int64, order_by::Symbol, 
    p::Parameters, maxsteps::Int64, option::Int64; hamming::Bool=false)
  sample_size_multiplier = 8
  sdf = DataFrame()
  sdf.goals = rand(MyInt(0):MyInt(size(df)[1]-1),sample_size_multiplier*sample_size)
  #sdf.goals = map(x->@sprintf("0x%x",x), rand(1:size(df)[1],sample_size_multiplier*sample_size))
  #println("sdf.goals: ",sdf.goals)
  #sdf[!,order_by] = [ df[ gg, order_by ] for gg in randgoallist(sample_size_multiplier*sample_size,p.numinputs,p.numoutputs) ] 
  sdf[!,order_by] = [df[ df.goal.==@sprintf("0x%x",gg),order_by] for gg in sdf.goals ]
  nsdf = sort(sdf,order(order_by))
  simple_goals = nsdf.goals[1:sample_size]
  complex_goals = nsdf.goals[((sample_size_multiplier-1)*sample_size+1):(sample_size_multiplier*sample_size)]
  steps_list = Int64[]
  steps_count = 0
  steps_squared_count = 0
  for i = 1:nreps
    if option == 1
      (src_g,dst_g) = choose_g_pair(simple_goals,simple_goals, sample_size, p, hamming=hamming)
    elseif option == 2
      (src_g,dst_g) = choose_g_pair(simple_goals,complex_goals, sample_size, p, hamming=hamming)
    elseif option == 3
      (src_g,dst_g) = choose_g_pair(complex_goals,simple_goals, sample_size, p, hamming=hamming)
    elseif option == 4
      (src_g,dst_g) = choose_g_pair(complex_goals,complex_goals, sample_size, p, hamming=hamming)
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

function choose_g_pair( goals1::Vector{MyInt}, goals2::Vector{MyInt}, sample_size::Int64, p::Parameters; hamming::Bool=false )
  if !hamming
    return (rand(goals1), rand(goals1))
  else 
    if rand() < 0.5
      choice1 = rand(goals1)
      cdf= DataFrame()
      cdf.goals2 = goals2
      cdf.chooseby = [ hamming_distance( choice1, gg, p.numinputs ) for gg in goals2 ] 
      sort!(cdf, [:chooseby])
      choice2 = rand(cdf.goals2[1:Int(ceil(sample_size/4))])
    else
      choice2 = rand(goals2)
      cdf= DataFrame()
      cdf.goals1 = goals1
      cdf.chooseby = [ hamming_distance( choice2, gg, p.numinputs ) for gg in goals1 ] 
      sort!(cdf, [:chooseby])
      choice1 = rand(cdf.goals1[1:Int(ceil(sample_size/4))])
    end
    return (choice1,choice2)
  end
end

# First, evolve a circuit that outputs goal s_g.  Then starting with this circuit,
#   evolve another circuit that outputs d_g.
function evolve_g_pair( s_g::Goal, d_g::Goal, p::Parameters, maxsteps::Int64 )
  funcs = default_funcs(p.numinputs)
  c = random_chromosome( p, funcs)
  (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) =
      mut_evolve_increase_numints(c, [s_g], funcs, maxsteps ) 
  println("src chromsome evolved in ",step," steps.")
  (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) =
      mut_evolve_increase_numints(c, [d_g], funcs, maxsteps ) 
  println("dst chromsome evolved in ",step," steps.")
  step
end

# Test the Wagner (2008) that genotypic robustness is negatively associated with genotypic evolvability.  
# Strongly confirmed by:
# julia> tge=test_g_evolvability_g_robustness( 100000, 4, 1, 9, 5 );
# julia> corspearman( tge[1], tge[2])
#  -0.8104812361089767
#=
function test_g_evolvability_g_robustness( nreps::Int64, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64 )
  p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
  funcs = default_funcs(numinputs)
  robust_vec = Float64[]
  all_outputs = Goal[]
  for i = 1:nreps
    c = random_chromosome(p,funcs)
    c_output = output_values(c)
    #println("c_output: ",c_output)
    outputs = mutate_all( c, funcs, output_outputs=true )
    #println("outputs: ",outputs)
    robust_outputs = filter( x->x==c_output, outputs )
    #println("robust_outputs: ",robust_outputs)
    evolvable_outputs = unique(outputs)
    all_outputs = unique( vcat( all_outputs, evolvable_outputs ) 
    #println("evolvable_outputs: ",evolvable_outputs)
    #println("push: ",length(robust_outputs)/length(outputs),"  ",length(evolvable_outputs)/length(outputs))
    push!( robust_vec, length(robust_outputs)/length(outputs))
    push!( evolvable_vec, length(evolvable_outputs)/length(outputs))
  end
  ( robust_vec, evolvable_vec )
end
=#

# Compute average genotypic robustness, average genotypic evolvability, and phenotypic evolvability for each of ngoals random goals.
# maxreps is attempted number of repetitions.  max_steps the the upper bound on the number of attempts including failures of mut_evolve().
function run_geno_robustness( ngoals::Int64, maxreps::Int64, numinputs::Int64, numoutputs::Int64, 
    numinteriors::IntRange, numlevelsback::Int64, max_steps::Int64, max_tries::Int64; csvfile = "" )
  if max_tries < maxreps
    error("max_tries should be greater than maxreps in geno_complexity.")
  end
  p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=10, numlevelsback=numlevelsback ) # Establish scope for p
  robust_vec = Float64[]
  evolvable_vec = Float64[]
  list_goals_params = Tuple[]
  for j = 1:ngoals
    goal = randgoal( numinputs, numoutputs)
    for numints in numinteriors 
      p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numints, numlevelsback=numlevelsback )
      println("goal: ",goal)
      push!(list_goals_params,(goal,p) )
    end
  end
  println("list goals_params: ",map(x->x[1],list_goals_params))
  result = pmap(x->geno_robustness( x[1], maxreps, p, max_steps, max_tries ), list_goals_params )
  goal_vec = [ rb[1] for rb in result ]
  robust_vec = [ rb[2] for rb in result ]
  geno_evolvable_vec = [ rb[3] for rb in result ]
  pheno_evolvable_vec = [ rb[4] for rb in result ]
  numints_vec = [ rb[5] for rb in result ] 
  ntries_vec = [ rb[6] for rb in result ]  
  nrepeats_vec = [ rb[7] for rb in result ]
  #println("(robust_vec, evolvable_vec): ",(robust_vec, evolvable_vec))
  robust_evo_df = DataFrame() 
  robust_evo_df.goal = goal_vec 
  robust_evo_df.robust = robust_vec
  robust_evo_df.geno_evolvable = geno_evolvable_vec
  robust_evo_df.pheno_evolvable = pheno_evolvable_vec
  robust_evo_df.numints = numints_vec
  robust_evo_df.ntries = ntries_vec
  robust_evo_df.nrepeats = nrepeats_vec
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )    
      #println(f,"# run time in minutes: ",ttime/60)
      println(f,"# funcs: ", Main.CGP.default_funcs(numinputs))
      println(f,"# ngoals: ",ngoals)
      println(f,"# maxreps: ",maxreps)
      println(f,"# max_tries: ",max_tries)
      println(f,"# numinputs: ",numinputs)
      println(f,"# numoutputs: ",numoutputs)
      println(f,"# numlevelsback: ",numlevelsback)
      println(f,"# max_steps: ",max_steps)
      CSV.write(f, robust_evo_df, append=true, writeheader=true )
    end
  end
  return robust_evo_df
end

# Compute average genotypic robustness, average genotypic evolvability, and phenotypic evolvability for goal.
# Return a row of the dataframe constructed in run_geno_robustness()
# maxreps is attempted number of repetitions.  max_steps the the upper bound on the number of attempts including failures of mut_evolve().
function geno_robustness( goal::Goal, maxreps::Int64, p::Parameters, max_steps::Int64, max_tries::Int64 )
  #p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
  funcs = default_funcs(p.numinputs)
  nrepeats = 0
  all_outputs_sum = 0
  robust_sum = 0
  g_evolvable_sum = 0
  all_unique_outputs = Goal[]
  i = 0
  while i < max_tries && nrepeats < maxreps
    i += 1
    #println("i: ",i,"  nrepeats: ",nrepeats)
    c = random_chromosome(p,funcs)
    c_output = output_values(c)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
      mut_evolve( c, [goal], funcs, max_steps )
    if step == max_steps
      println("mut evolve failed for goal: ",goal)
      continue
    end
    c_output = output_values(c)
    @assert c_output == goal
    #println("output values: ",c_output)
    outputs = mutate_all( c, funcs, output_outputs=true )
    all_outputs_sum += length(outputs)
    robust_outputs = filter( x->x==c_output, outputs )
    robust_sum += length(robust_outputs)
    p_evolvable_outputs = unique(outputs)
    g_evolvable_sum += length(p_evolvable_outputs)
    all_unique_outputs = unique(vcat(all_unique_outputs,p_evolvable_outputs))
    nrepeats += 1
    #println((all_outputs_sum,robust_sum,length(all_unique_outputs)))
  end  
  ntries = i
  if nrepeats > 0

  else
    return ( goal, 0.0, 0.0, 0.0, p.numinteriors, ntries, nrepeats )
  end
end

# For each of the goals in goallist, evolves maxreps chromosomes whose output is that goal.
# Computes properties of these goals and Chromosomes/LinCircuits.. 
# Arguments:
#   maxreps:  The number of evolutions attempted for each goal
#   iter_maxreps:  The number of evolutions attemped in each pmap-parallel iteration
#   max_steps:  The maximum number of steps for each run of neutral_evolution.
#   max_tries:  The maximum number of evolutions tried in each pmap-parallel iteration 
#   consolidate:  if true, consolidate multiple rows with the same goal
# Some goals are much more difficult to evolve than others.  To spread the evolutions of difficult
#   goals over more pmap-parallel iterations, set iter_maxreps to be small relative to maxreps
# The max_tries parameter determines how many attempts are made to evolve difficult goals.
# Returns a dataframe with num_iterations=Int(ceil(maxreps/iter_maxreps)) rows per goal.
function run_geno_complexity( goallist::GoalList, maxreps::Int64, iter_maxreps::Int64, p::Parameters,
      max_steps::Int64, max_tries::Int64, maxsteps_recover::Int64=0, maxtrials_recover::Int64=0, maxtries_recover::Int64=0; 
      use_lincircuit::Bool=false, consolidate::Bool=true, csvfile::String = "" )
  #p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=10, numlevelsback=numlevelsback ) # Establish scope for p
  geno_complexity_df = DataFrame() 
  geno_complexity_df.goal = Vector{MyInt}[]
  geno_complexity_df.numinputs = Int64[]
  geno_complexity_df.numoutputs = Int64[]
  if use_lincircuit
    geno_complexity_df.numinstructs = Int64[]
    geno_complexity_df.numregisters = Int64[]
  else
    geno_complexity_df.numints = Int64[]
    geno_complexity_df.numlevsback = Int64[]
  end
  geno_complexity_df.maxsteps = Int64[]
  geno_complexity_df.ntries = Int64[]
  geno_complexity_df.nsuccesses = Float64[]
  geno_complexity_df.avg_steps = Float64[]
  geno_complexity_df.log_avg_steps = Float64[]
  geno_complexity_df.robustness = Float64[]
  geno_complexity_df.evo_count = Int64[]
  geno_complexity_df.unique_goals = GoalList[]
  geno_complexity_df.nactive = Float64[]
  geno_complexity_df.complexity = Float64[]
  geno_complexity_df.degeneracy = Float64[]
  geno_complexity_df.sumsteps = Float64[]
  geno_complexity_df.sumtries = Float64[]
  #geno_complexity_df.complexQ95 = Float64[]
  #geno_complexity_df.complexQ99 = Float64[]
  #geno_complexity_df.epi2= Float64[]
  #geno_complexity_df.epi3= Float64[]
  #geno_complexity_df.epi4= Float64[]
  #geno_complexity_df.epi_total = Float64[]
  #geno_complexity_df.f_mutrobust= Float64[]
  list_goals = Goal[]
  num_iterations = Int(ceil(maxreps/iter_maxreps))  # The number of iterations for each goal
  iter_maxtries = Int(ceil(max_tries/num_iterations))
  println("num_iterations: ",num_iterations,"  iter_maxtries: ",iter_maxtries,"  iter_maxreps: ",iter_maxreps)
  if iter_maxtries < iter_maxreps
    error("iter_maxtries should be greater than iter_maxreps in run_geno_complexity.")
  end
  for g in goallist
    for i = 1:num_iterations
      #println("i: ",i,"  goal: ",g)
      push!(list_goals,g)
    end
  end
  funcs = default_funcs(p.numinputs)
  c = random_chromosome(p,funcs)
  sample_size = length(mutate_all(c, funcs, output_outputs=true))*iter_maxreps 
  num_goals = 2^2^p.numinputs
  #println("sample size: ",sample_size,"  num_goals: ",num_goals)
  #println("list_goals: ",list_goals)
  result = pmap(g->geno_complexity( g, iter_maxreps, p, max_steps, iter_maxtries,
      maxsteps_recover, maxtrials_recover, iter_maxtries, use_lincircuit=use_lincircuit ), list_goals)
  #result = map(g->geno_complexity( g, iter_maxreps, p, max_steps, iter_maxtries,
  #    maxsteps_recover, maxtrials_recover, maxtries_recover, use_lincircuit=use_lincircuit ), list_goals)
  #println("after pmap:  size(geno_complexity_df): ",size(geno_complexity_df))
  for res in result
    push!(geno_complexity_df,res)
  end
  geno_complexity_df.evo_count = zeros(Int64,size(geno_complexity_df)[1] )
  #geno_complexity_df.ratio = zeros(Float64,size(geno_complexity_df)[1] )
  #geno_complexity_df.estimate = zeros(Float64,size(geno_complexity_df)[1] )
  j = 1
  for g in goallist
    println("goal ",g)
    all_unique_goals = Goal[]
    prev_evo_count = 0
    for i = 1:num_iterations
      all_unique_goals = unique( vcat( all_unique_goals, geno_complexity_df[j,:unique_goals] ))
      evo_count = length(all_unique_goals)
      #println("goal: ",geno_complexity_df[j,:goal],"  j: ",j,"  i: ",i,"  length(all_unique_goals): ",length(all_unique_goals))
      new_count = evo_count - prev_evo_count
      geno_complexity_df[j,:evo_count] = evo_count
      prev_evo_count = evo_count
      j += 1
    end
  end
  select!(geno_complexity_df,DataFrames.Not(:unique_goals))    # Remove :unique_goals from gcdf
  if consolidate   # Consolidate multiple rows with the same goal.  See Analyze.jl for code.
    geno_complexity_df = consolidate_dataframe( geno_complexity_df )
  end
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )    
      #println(f,"# run time in minutes: ",ttime/60)
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      println(f,"# use_lincircuit: ",use_lincircuit)
      println(f,"# maxreps: ",maxreps)
      println(f,"# iter_maxreps: ",iter_maxreps)
      println(f,"# numinputs: ",p.numinputs)
      println(f,"# numoutputs: ",p.numoutputs)
      if use_lincircuit
        println(f,"# numinstructions: ",p.numinteriors) 
        println(f,"# numregisters: ",p.numlevelsback)
      else
        println(f,"# numinteriors: ",p.numinteriors)
        println(f,"# numlevelsback: ",p.numlevelsback)
      end
      println(f,"# max_steps: ",max_steps)
      CSV.write(f, geno_complexity_df, append=true, writeheader=true )
    end
  end
  return geno_complexity_df
end

# For the given goal, evolves iter_maxreps chromosomes that output that goal.
# Returns a tuple which is pushed as a row onto the dataframe constructed in run_geno_complexity().
# This can be run in parallel for one goal because each run is evolving iter_maxreps circuits that compute goal.
# all_unique_outputs is returned as one of the fields in the returned dataframe.
# Then run_geno_complexity combines the all_unique_outputs from the parallell runs
function geno_complexity( goal::Goal, iter_maxreps::Int64, p::Parameters,  maxsteps::Int64, max_tries::Int64,
      maxsteps_recover::Int64, maxtrials_recover::Int64, maxtries_recover::Int64; use_lincircuit::Bool=false )
  #println("geno_complexity: goal: ",goal)
  funcs = default_funcs(p.numinputs)
  #W = Walsh(2^p.numinputs)
  all_outputs_sum = 0
  robust_sum = 0
  all_unique_outputs = Goal[]
  nactive_list = Int64[]
  complexity_list = Float64[]
  degeneracy_list = Float64[]
  sumsteps_list = Int64[]
  sumtries_list = Int64[]
  #frenken_mi_list = Float64[]
  if p.numoutputs == 0
    epi2 = k_bit_epistasis(W,2,goal[1])
    epi3 = k_bit_epistasis(W,3,goal[1])
    epi4 = k_bit_epistasis(W,4,goal[1])
    epi_total = total_epistasis(W,goal[1])
  end
  #println("After epi")
  sum_steps = 0.0
  sum_steps_per_iteration = 0  # for one iteration of the while loop
  numsuccesses = 0  # The number of successful evolutions
  i = 0
  # Do up to max_tries evolutions attempting to do iter_maxreps successful evolutions
  while i < max_tries && numsuccesses < iter_maxreps 
    i += 1
    c = use_lincircuit ? rand_lcircuit(p,funcs) : random_chromosome(p,funcs)
    #println("i: ",i,"  funcs: ",funcs)
    #print_circuit(c,funcs)
    c_output = output_values(c)
    #(c,steps,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [goal], funcs, maxsteps )
    (c,steps) = neutral_evolution( c, funcs, goal, maxsteps )
    if steps == maxsteps
      #println("mut evolve failed for goal: ",goal)
      #println("neutral evolution failed for goal: ",goal)
      sum_steps_per_iteration += steps
      continue
    end        
    c_output = output_values(c)
    #println("goal: ",goal,"  i: ",i,"  numsuccesses: ",numsuccesses,"  c_output: ",c_output)
    @assert sort(c_output) == sort(goal)
    #println("output values: ",c_output,"  complexity5(c): ",complexity5(c))
    outputs = mutate_all( c, funcs, output_outputs=true )
    all_outputs_sum += length(outputs)
    robust_outputs = filter( x->x==c_output, outputs )
    robust_sum += length(robust_outputs)
    evolvable_outputs = unique(outputs)
    all_unique_outputs = unique(vcat(all_unique_outputs,evolvable_outputs))
    push!( nactive_list, number_active( c ))
    complexity = use_lincircuit ? lincomplexity( c, funcs ) : complexity5(c)
    push!( complexity_list, complexity )
    push!( degeneracy_list, degeneracy( c ))
    (sumsteps,sumtries) = recover_phenotype( c, maxsteps_recover, maxtrials_recover, maxtries_recover )
    push!( sumsteps_list, sumsteps )
    push!( sumtries_list, sumtries )
    #push!( frenken_mi_list, fmi_chrome( c ))
    sum_steps_per_iteration += steps
    sum_steps += sum_steps_per_iteration
    sum_steps_per_iteration = 0
    numsuccesses += 1
  end  
  ntries = i
  #println("ntries: ",ntries,"  numsuccesses: ",numsuccesses,"  all_outputs_sum: ",all_outputs_sum,"  len all_unique_outputs: ",length(all_unique_outputs))
  println("num_successes: ",numsuccesses)
  # Return evolvability count for testing  10/13
  #length(all_unique_outputs)
  # Temporarily comment out for testing evolvability  10/13
  if numsuccesses > 0
    return ( 
      goal, 
      p.numinputs,
      p.numoutputs,
      p.numinteriors,
      p.numlevelsback,
      maxsteps,
      ntries,
      numsuccesses,
      sum_steps/numsuccesses,   # Note:  Only counts steps for successful runs.
      sum_steps != 0.0 ? log10(sum_steps/numsuccesses) : 0.0,
      robust_sum/all_outputs_sum, 
      0,   # evo_count, value filled in later
      #0.0,  # ratio, value filled in later 
      #0.0,  # estimate, value filled in later 
      all_unique_outputs,
      sum( nactive_list )/iter_maxreps,
      sum( complexity_list )/iter_maxreps,
      sum( degeneracy_list )/iter_maxreps,
      sum( sumsteps_list )/iter_maxreps,
      sum( sumtries_list )/iter_maxreps,
      #quantile(complexity_list,0.95),
      #quantile(complexity_list,0.99),
      #p.numoutputs==0 ? epi2 : 0.0,
      #p.numoutputs==0 ? epi3 : 0.0,
      #p.numoutputs==0 ? epi4 : 0.0,
      #p.numoutputs==0 ? epi_total : 0.0,
      #sum(frenken_mi_list)/iter_maxreps
    )
  else  # Evolution always failed
    return ( 
      goal,
      p.numinputs,
      p.numoutputs,
      p.numinteriors,
      p.numlevelsback,
      maxsteps,
      ntries,
      0,   # numsuccesses
      0,   # sum_steps
      0,   # log_sum_steps
      0.0, # robust sum
      0,   # evo_count
      #0.0, # ratio
      #0.0, # estimate
      Goal[],
      sum( nactive_list )/iter_maxreps,
      sum( complexity_list )/iter_maxreps,
      0.0,  # degeneracy
      0.0,  #sumsteps_list
      0.0,  #sumtries_list
      #0.0,  # quantile(complexity_list,0.95)
      #0.0,
      #0.0,
      #0.0,
      #0.0,
      #0.0
    )
  end
end

function add_frequencies_to_dataframe( gdf:: DataFrame, counts_field::Symbol, count_csv_file::String="../data/counts/count_out_4x1_all_ints_10_10.csv" )
  cdf = read_dataframe(count_csv_file)
  # The following is what works when the :goal field of gdf is of the format "UInt16[0x4fd1]"  where UInt16 is MyInt
  #counts = [cdf[cdf.goals.==@sprintf("0x%x",eval(Meta.parse(gdf.goal[i]))[1]),counts_field][1] for i = 1:size(gdf)[1]]
  # The following is what works when the :goal field of gdf is of the format [0x4fd1] where UInt16 is MyInt
  counts = [cdf[cdf.goals.==@sprintf("0x%x",gdf.goal[i][1]),counts_field][1] for i = 1:size(gdf)[1]]
  gdf.counts = counts
  gdf
end

# Create and save a scatter plot using Plots package
# For examples, see notes/diary10_11.txt
function scatter_plot( gcdf::DataFrame, y_var::Symbol, x_var::Symbol, circuit_type::String, numints::Int64, numlevsback::Int64 )
  Plots.scatter( gcdf[!,x_var], gcdf[!,y_var], title="$y_var vs $x_var $circuit_type $numints ints $numlevsback levsback", ylabel=y_var, xlabel=x_var, label="")
  fname = "$y_var@vs@$x_var@$circuit_type@$numints@ints@$numlevsback@levsback"  # "@" will be changed to "_"
  Plots.savefig( replace(fname, "@"=>"_") )
  # Redo plot so that it shows up interactively
  Plots.scatter( gcdf[!,x_var], gcdf[!,y_var], title="$y_var vs $x_var $circuit_type $numints ints $numlevsback levsback", ylabel=y_var, xlabel=x_var, label="")
end

# See diary9_26.txt for motivation and sample results
# TODO:  filter parent chromosomes for high or low complexity
function parent_child_complexity( p::Parameters, nreps::Int64 )
  nrepeats = 3
  funcs = default_funcs(p.numinputs)
  complexity_pair_list = Tuple{Float64,Float64}[]
  for i = 1:nreps
    pc = random_chromosome(p,funcs)
    p_complex = complexity5( pc )
    cc = mutate_chromosome!( deepcopy(pc), funcs )[1]
    c_complex = complexity5( cc )
    #println("pair: ",(p_complex, c_complex ))
    push!( complexity_pair_list, (p_complex, c_complex ) )
  end
  complexity_pair_list
end

rand_norm( mean::Float64, std::Float64, nreps::Int64 ) = std.*randn(nreps).+mean
rand_norm(nreps::Int64) = rand_norm( 3.12727, 1.164145, nreps )

# population is a vector of chromsomes
# Returns the pair (evolvability,robustness)
# The evolvability of the population is the length of the set of unique phenotypes in the union 
#  of mutational neighborhoods of the genotypes of the population 
#
function pop_evolvability_robustness( population::Vector{Chromosome} )
  p = population[1].params
  funcs = default_funcs(p.numinputs)
  robust_sum = 0.0
  count_neighbors = 0
  ph_set = Set(Vector{Vector{MyInt}}())
  for ch in population
    neighbors = mutate_all( ch, funcs, output_outputs=true )
    count_neighbors += length(neighbors)
    robustness = length(filter(x->x==output_values(ch),neighbors))/length(neighbors)
    robust_sum += robustness
    for neighbor in neighbors 
      push!( ph_set, neighbor )
    end
  end
  (length(ph_set)/count_neighbors,robust_sum/length(population))
end
#=
# Uses map-reduce to union the set results of multiple runs of pheno_set_rand_neutral_walk() for a single phenotype ph.
# unions of the results of multiple walks
function pheno_set_rand_neutral_walks( p::Parameters, funcs::Vector{Func}, ph::Goal, nwalks::Int64, walk_length::Int64,  max_tries::Int64, max_steps::Int64;
    use_lincircuit::Bool=false, csvfile::String="" )
  # One list of common phenotypes included here for reference.
  #common_list = [ 0x0000, 0x0011, 0x0022, 0x0033, 0x0044, 0x0055, 0x0066, 0x0077, 0x0088, 0x0099, 0x00aa, 0x00bb, 0x00cc, 0x00dd, 0x00ee, 0x00ff ]
  pheno_sets = pmap( _->pheno_set_rand_neutral_walk( p, funcs, ph, walk_length, max_tries, max_steps, use_lincircuit=use_lincircuit ), collect(1:nwalks ) )
  pheno_set = reduce( union, pheno_sets )
  df = DataFrame()
  df.pheno_list =  [ s for s in pheno_set ] 
  if length(csvfile) > 0
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  open( csvfile, "w" ) do f
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    #println(f,"# run time in minutes: ",(ptime+ntime)/60)
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ", Main.CGP.default_funcs(p))
    println(f,"# phenotype: ", ph)
    println(f,"# nwalks: ",nwalks)
    println(f,"# walk_length: ",walk_length)
    println(f,"# max_tries: ",max_tries)
    println(f,"# max_steps: ",max_steps)
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

# returns the set of phenotypes encountered in a random neutral walk starting from a circuit evolved to map to phenotype ph.
function pheno_set_rand_neutral_walk( p::Parameters, funcs::Vector{Func}, ph::Goal, walk_length::Int64,  max_tries::Int64, max_steps::Int64;
    use_lincircuit::Bool=false )
  (nc,steps) = pheno_evolve( p, funcs, ph::Goal, max_tries, max_steps; use_lincircuit=use_lincircuit ) 
  if steps == max_steps
    println("pheno_evolve failed to evolve a circuit to map to phenotype: ", ph )
    return Set(MyInt[])
  end
  pheno_set_rand_neutral_walk( nc, funcs, walk_length, max_tries, max_steps )
end

# returns the set of phenotypes encountered in a random neutral walk starting from a circuit evolved to map to phenotype ph.
function pheno_set_rand_neutral_walk( c::Circuit, funcs::Vector{Func}, walk_length::Int64,  max_tries::Int64, max_steps::Int64 )
  use_lincircuit = typeof(c) == LinCircuit
  println("use_lincircuit: ",use_lincircuit)
  p = c.params
  @assert p.numoutputs == 1
  goal = output_values(c)[1]
  pheno_set = Set(MyInt[])  # Assumes 1 output
  for i = 1:walk_length
    j = 1
    while j <= max_steps   # also terminated by a break statement
      sav_c = deepcopy(c)
      use_lincircuit ? mutate_circuit!( c, funcs ) : mutate_chromosome!( c, funcs )
      #println("c: ",c)
      output = output_values(c)[1]
      if output == goal
        #println("successful step for i= ",i,"  outputs: ",outputs)
        all_neighbor_phenos = mutate_all( c, funcs )   # polymorphic call
        #println("j: ",j,"  all_neighbor_phenos: ",all_neighbor_phenos)
        pheno_list = map( x->x[1], all_neighbor_phenos )
        #println("j: ",j,"  pheno_list: ",pheno_list)
        union!(pheno_set, Set(pheno_list) )
        break
      end
      c = sav_c
      j += 1
    end
    if step == max_steps
      error("Failed to find neutral mutation")
    end
  end
  pheno_set
end
=#
