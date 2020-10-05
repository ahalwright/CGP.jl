# Define phenotypic evolvability where a phenotype is a goal or chromosome output.
# Note that genotypic evolvability is computed by calling mutate_all( c, funcs, robust_only=true )[2]
# Use Wagner's (2008) method of evolving nchromess chromsomes whose output is the goal (where Wagner uses nchromes=100).

using DataFrames
using Dates
using CSV
using Distributed
using HypothesisTests
using Printf

export evolvability, run_evolvability, evo_result, test_evo, evo_result_type, run_evolve_g_pairs 
export run_geno_robustness, geno_robustness
#export run_geno_robustness0, geno_robustness0
export run_geno_complexity, geno_complexity
export parent_child_complexity, rand_norm
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

function run_geno_robustness( ngoals::Int64, maxreps::Int64, numinputs::Int64, numoutputs::Int64, 
    numinteriors::IntRange, numlevelsback::Int64, max_steps::Int64, max_tries::Int64; csvfile = "" )
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
  evolvable_vec = [ rb[3] for rb in result ]
  numints_vec = [ rb[4] for rb in result ] 
  ntries_vec = [ rb[5] for rb in result ]  
  nrepeats_vec = [ rb[6] for rb in result ]
  #println("(robust_vec, evolvable_vec): ",(robust_vec, evolvable_vec))
  robust_evo_df = DataFrame() 
  robust_evo_df.goal = goal_vec 
  robust_evo_df.robust = robust_vec
  robust_evo_df.evolvable = evolvable_vec
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

function geno_robustness( goal::Goal, maxreps::Int64, p::Parameters, max_steps::Int64, max_tries::Int64 )
  #p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
  funcs = default_funcs(p.numinputs)
  nrepeats = 0
  all_outputs_sum = 0
  robust_sum = 0
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
    evolvable_outputs = unique(outputs)
    all_unique_outputs = vcat(all_unique_outputs,evolvable_outputs)
    nrepeats += 1
  end  
  ntries = i
  if nrepeats > 0
    return ( goal, robust_sum/all_outputs_sum, length(all_unique_outputs)/all_outputs_sum, p.numinteriors, ntries, nrepeats )
  else
    return ( goal, 0.0, 0.0, p.numinteriors, ntries, nrepeats )
  end
end

# For each of ngoals randomly chosen goals, evolves nreps chromosomes whose output is that goal.
# Computes properties of these goals and chromosomes. 
# Returns a dataframe with one row per goal.
function run_geno_complexity( ngoals::Int64, nreps::Int64, numinputs::Int64, numoutputs::Int64, 
    numinteriors::IntRange, numlevelsback::Int64, max_steps::Int64; csvfile = "" )
  robust_vec = Float64[]
  evolvable_vec = Float64[]
  list_goals_params = Tuple[]
  for j = 1:ngoals
    goal = randgoal( numinputs, numoutputs)
    for numints in numinteriors 
      #println("goal: ",goal)
      p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
      push!(list_goals_params,(goal,p) )
    end
  end
  #result = pmap(x->geno_complexity( x[1], nreps, x[2], max_steps ), list_goals_params )
  result = map(x->geno_complexity( x[1], nreps, x[2], max_steps ), list_goals_params )
  #println("(robust_vec, evolvable_vec): ",(robust_vec, evolvable_vec))
  geno_complexity_df = DataFrame() 
  geno_complexity_df.goal = MyInt[]
  geno_complexity_df.numinputs = Int64[]
  geno_complexity_df.numoutputs = Int64[]
  geno_complexity_df.numints = Float64[]
  geno_complexity_df.numlevsback = Float64[]
  geno_complexity_df.maxsteps = Int64[]
  geno_complexity_df.robustness = Float64[]
  geno_complexity_df.evolvability = Float64[]
  geno_complexity_df.nactive = Float64[]
  geno_complexity_df.complexity = Float64[]
  geno_complexity_df.epi2= Float64[]
  geno_complexity_df.epi3= Float64[]
  geno_complexity_df.epi4= Float64[]
  geno_complexity_df.epi_total = Float64[]
  geno_complexity_df.f_mutrobust= Float64[]
  for res in result
    push!(geno_complexity_df,res)
  end
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )    
      #println(f,"# run time in minutes: ",ttime/60)
      println(f,"# funcs: ", Main.CGP.default_funcs(numinputs))
      println(f,"# ngoals: ",ngoals)
      println(f,"# nreps: ",nreps)
      println(f,"# numinputs: ",numinputs)
      println(f,"# numoutputs: ",numoutputs)
      println(f,"# numinteriors: ",numinteriors)
      println(f,"# numlevelsback: ",numlevelsback)
      println(f,"# max_steps: ",max_steps)
      CSV.write(f, geno_complexity_df, append=true, writeheader=true )
    end
  end
  return geno_complexity_df
end

# For the given goal, evolves nreps chromosomes that output that goal.
# Returns a tuple which is pushed as a row onto the dataframe constructed in run_geno_complexity().
function geno_complexity( goal::Goal, nreps::Int64, p::Parameters,  maxsteps::Int64 )
  #p = Parameters( numinputs = numinputs, numoutputs = numoutputs, numinteriors = numinteriors, numlevelsback = numlevelsback )
  funcs = default_funcs(p.numinputs)
  W = Walsh(2^p.numinputs)
  all_outputs_sum = 0
  robust_sum = 0
  all_unique_outputs = Goal[]
  new_numints = p.numinteriors
  new_numlevback = p.numlevelsback
  new_numints_list = zeros(Int64, nreps )
  new_levsback_list = zeros(Int64, nreps )
  nactive_list = zeros(Int64, nreps )
  complexity_list = zeros(Float64, nreps )
  frenken_mi_list = zeros(Float64, nreps ) 
  epi2 = k_bit_epistasis(W,2,goal[1])
  epi3 = k_bit_epistasis(W,3,goal[1])
  epi4 = k_bit_epistasis(W,4,goal[1])
  epi_total = total_epistasis(W,goal[1])
  for i = 1:nreps
    #new_numints = p.numinteriors
    c = random_chromosome(p,funcs)
    c_output = output_values(c)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback) =
      mut_evolve_increase_numints( c, [goal], funcs, maxsteps )
    c_output = output_values(c)
    @assert c_output == goal
    #println("output values: ",c_output)
    outputs = mutate_all( c, funcs, output_outputs=true )
    all_outputs_sum += length(outputs)
    robust_outputs = filter( x->x==c_output, outputs )
    robust_sum += length(robust_outputs)
    evolvable_outputs = unique(outputs)
    all_unique_outputs = vcat(all_unique_outputs,evolvable_outputs)
    nactive_list[i] = number_active( c )
    complexity_list[i] = complexity5( c )
    new_numints_list[i] = new_numints
    new_levsback_list[i] = new_levsback
    frenken_mi_list[i] = fmi_chrome( c )
  end  
  return ( 
    goal[1], 
    p.numinputs,
    p.numoutputs,
    sum( new_numints_list )/nreps,
    sum( new_levsback_list )/nreps,
    maxsteps,
    robust_sum/all_outputs_sum, 
    length(all_unique_outputs)/all_outputs_sum,
    sum( nactive_list )/nreps,
    sum( complexity_list )/nreps,
    epi2,
    epi3,
    epi4,
    epi_total,
    sum(frenken_mi_list)/nreps
  )
end

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

function test_evo()
  g = [0x0332]
  numinputs=4; numoutputs=1; numinteriors=9; numlevelsback=5;
  setup_funcs(numinputs)
  g =[0x03f3]
  p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numinteriors, numlevelsback=numlevelsback )
  max_steps = 200000
  nchromes = 50
  nrepeats = 3
  evo_result = Main.CGP.evo_result
  er = evo_result( g, nchromes, max_steps, p, nrepeats, 0, 0, 0 )
  funcs=default_funcs(p.numinputs)
  #evolvability( er, funcs )
  #evolvability( er, funcs, intermediate_gens=[1,4,8] )
  run_evolvability( 2, [g], funcs, er.nchromes, er.maxsteps, er.numinputs, er.numints, er.levelsback, er.nrepeats, intermediate_gens=[10,20,30,40,50] )
  #run_evolvability( 20, [g], funcs, 5:5:20, er.maxsteps, er.numinputs, er.numints, er.levelsback, er.nrepeats, "test.csv" )
end

function test_gen_evo()
  df = read_dataframe("../data/consolidate/geno_pheno_raman_df_all_9_13.csv");
  p = Parameters( numinputs=4, numinteriors=10, numlevelsback=5, numoutputs=1 ); maxsteps = 150000
  sample_size = 100
  nreps = 100
  res = sample_g_pairs( df, 4, 10, :complex, p, maxsteps, 2 ) 
  @time result = run_sample_g_pairs( df, nreps, sample_size, :complex, p, maxsteps )
end


