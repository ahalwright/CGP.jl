# Define phenotypic evolvability where a phenotype is a goal or chromosome output.
# Note that genotypic evolvability is computed by calling mutate_all( c, funcs, robust_only=true )[2]
# Use Wagner's (2008) method of evolving nchromess chromsomes whose output is the goal (where Wagner uses nchromes=100).

using DataFrames
using Dates
using CSV
using Distributed

export evolvability, run_evolvability, test_evo

mutable struct evo_result_type
  goal::Goal
  nchromes::Int64
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  maxsteps::Int64
  evolvability::Int64
end

function evo_result( goal::Goal, nchromes, max_steps::Int64, p::Parameters, evolvability::Int64 )
  evo_result_type(
    goal,
    nchromes,
    p.numinputs,
    p.numoutputs,
    p.numinteriors,
    p.numlevelsback,  
    max_steps,
    evolvability
  )
end

function evo_result( goal::Goal, nchromes, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64, 
      max_steps::Int64, evolvability::Int64 )
  evo_result_type(
    goal,
    nchromes,
    numinputs,
    numoutputs,
    numinteriors,
    numlevelsback,  
    max_steps,
    evolvability
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
    er.evolvability
  )
end

function evolvability( g::Goal, funcs::Vector{Func}, nchromes::Int64, maxsteps::Int64, 
      numinputs::Int64, numinteriors::Int64, numlevelsback::Int64 )
  p = Parameters( numinputs=numinputs, numoutputs=length(g), numinteriors=numinteriors, numlevelsback=numlevelsback )
  evolvability( g, funcs, nchromes,  maxsteps, p )
end

function evolvability( er::evo_result_type, funcs::Vector{Func} )
  repeat_limit = 10
  println("er: ",er)
  p = Parameters( numinputs=er.numinputs, numoutputs=er.numoutputs, numinteriors=er.numints, numlevelsback=er.levelsback )
  result = Goal[]
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
    #print_build_chromosome( c )R
    goal_list = mutate_all( c, funcs, output_outputs=true )
    #println("goal_list: ",goal_list)
    goal_set = unique(goal_list)
    result = unique(vcat(result,goal_set))
  end
  println("len(result): ",length(result))
  er.evolvability = length(result)
  return er
end

function run_evolvability( nreps::Int64, gl::GoalList, funcs::Vector{Func}, nchromes::IntRange,  maxsteps::IntRange, 
      numinputs::Int64, numinteriors::IntRange, numlevelsback::IntRange )
  evo_result_list = evo_result_type[]
  df = DataFrame()
  df.goal=Vector{MyInt}[]
  df.nchromes=Int64[]
  df.numinputs=Int64[]
  df.numints=Int64[]
  df.levsback=Int64[]
  df.maxsteps=Int64[]
  df.evolvability=Int64[]
  for g in gl
    for nch = nchromes
      for nints = numinteriors
        for levsback = numlevelsback
          for mxsteps = maxsteps
            for rep = 1:nreps
              #er = evo_result( g, nch, mxsteps, p, 0 )
              er = evo_result( g, nch, numinputs, length(g), nints, levsback, mxsteps, 0 )
              #er = evolvability( er, funcs )
              #push!(df,(g,nch,numinputs,nints,levsback,mxsteps,evob))
              push!(evo_result_list, er )
            end
          end
        end
      end
    end
  end
  new_evo_result_list = pmap( x->evolvability(x,funcs), evo_result_list )
  #new_evo_result_list = map( x->evolvability(x,funcs), evo_result_list )
  for er in new_evo_result_list
    new_row = evo_result_to_tuple(er)
    push!(df, new_row)
  end
  df
end

function run_evolvability( nreps::Int64, gl::GoalList, funcs::Vector{Func}, nchromes::IntRange,  maxsteps::IntRange, 
      numinputs::Int64, numinteriors::IntRange, numlevelsback::IntRange, csvfile::String )
  #result = @timed run_evolvability( nreps, gl, funcs, nchromes,  maxsteps, numinputs, numinteriors, numlevelsback ) 
  (df,ttime) = @timed run_evolvability( nreps, gl, funcs, nchromes,  maxsteps, numinputs, numinteriors, numlevelsback ) 
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

function test_evo()
  g = [0x0332]
  numinputs=4; numoutputs=1; numinteriors=9; numlevelsback=4;
  setup_funcs(numinputs)
  g =[0x03f3]
  p = Parameters( numinputs=numinputs, numoutputs=length(g), numinteriors=numinteriors, numlevelsback=numlevelsback )
  nchromes = 3
  max_steps = 10000
  er = evo_result( g, nchromes, max_steps, p, 0 )
  funcs=default_funcs(p.numinputs)
  evolvability( er, funcs )
  run_evolvability( 20, [g], funcs, 5:5:20, er.maxsteps, er.numinputs, er.numints, er.levelsback, "test.csv" )
end
