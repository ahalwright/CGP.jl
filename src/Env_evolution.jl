# Implement experiments to measure the selection for env. robustness evolution as a function of:
#  1.  Funcs
#  2.  numinteriors
#  3.  numinputs #  4.  numoutputs
#  5.  max_steps
# Uses randomly generated goallist
# Keep track of number of "worse" and "same" updates of Chromosome since its last improvement
export run_env_evolution, env_result_to_tuple
using DataFrames 
using CSV
using Statistics
using Distributions
using Distributed
env_result_type = Main.CGP.env_result_type
#=
iterations = 4
numinputs = 2:2
numoutputs = 2:2
nodearity = 2
numinteriors = 6:6
numlevelsback = 6:6
ngoals = 4:4
goallistlength=8:8
levelsback=6:6
maxsteps = 400:400
run_mut_evolution( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback, "testdata.csv" )
=#
#context = construct_contexts(numinputs)[numinputs]
#p = Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
#iterations = 20
#print_parameters(p)

function run_env_evolution( numiterations::Int64, numinputs::IntRange, numoutputs::IntRange, 
    numinteriors::IntRange, goallistlength::IntRange, maxsteps::IntRange,
    levelsback::IntRange, gl_repetitions::IntRange, num_flip_bits::IntRange, avgfitness::IntRange, 
    perm_heuristic::IntRange, perturb_goal_range::IntRange, fit_limit_list::Vector{Float64}, csvfile::String; 
    base::Float64=2.0, active_only::Bool=false )
  hamming_sel = true
  maxints_for_degen = 20
  nodearity = 2
  perturb_goal=true
  env_result_list = env_result_type[]
  df = DataFrame() 
  df.numinputs=Int64[]
  df.numoutputs=Int64[]
  df.numints=Int64[]
  df.levelsback=Int64[]
  df.ngoals=Int64[]
  #df.hamming_sel=Bool[]
  #df.robust_sel=Bool[]
  #df.active_only=Bool[]
  df.maxsteps=Int64[]
  df.gl_reps=Int64[]
  df.nflipbits=Int64[]
  df.perturb_goal=Bool[]
  df.avgfitness=Bool[]
  df.fitlimit=Float64[]
  df.steps=Int64[]
  df.avgfit=Float64[]
  #df.same=Int64[]
  #df.worse=Int64[]
  #df.better=Int64[]
  df.nactive=Int64[]
  df.redundancy=Float64[]
  df.complexity=Float64[]
  df.degeneracy=Float64[]
  df.sdegeneracy=Float64[]
  #println("size(df): ",size(df))
  for num_inputs = numinputs
    for num_outputs = numoutputs
      funcs = default_funcs(num_inputs) 
      for num_interiors = numinteriors
        for num_goals = goallistlength
          println("numinputs: ",num_inputs,"  numoutputs: ",num_outputs,"  numints: ",num_interiors,"  numgoals: ",num_goals)
          for levsback = levelsback
            for gl_reps = gl_repetitions 
              for avgfit = avgfitness
                for perturb_goal = perturb_goal_range
                  for nflipbits = num_flip_bits
                    #println("perturb_goal: ",perturb_goal)
                    for max_steps = maxsteps
                      for fit_limit in fit_limit_list
                        #println("iters fit_limit: ",fit_limit)
                        for _ = 1:numiterations
                          p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, levsback )
                          rr = env_result( p, num_goals, hamming_sel, active_only, max_steps, gl_reps, nflipbits, 
                              perm_heuristic, perturb_goal, avgfit, fit_limit )
                          push!(env_result_list,rr)
                          #run_mut_evolve!(rr)
                          #new_row = env_result_to_tuple(rr)
                          #Base.push!( df, new_row )
                        end
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
  new_env_result_list = pmap(r->run_env_evolve!(r,maxints_for_degen=maxints_for_degen,base=base),env_result_list)
  #new_env_result_list = map(r->run_env_evolve!(r,maxints_for_degen=maxints_for_degen,base=base),env_result_list)
  for r = new_env_result_list
     new_row = env_result_to_tuple(r)
     #println("length(new_row): ",length(new_row),"  size(df): ",size(df)) 
     Base.push!( df, new_row )
  end
  println(default_funcs(2))
  open( csvfile, "w" ) do f
    println(f,"funcs: ", Main.CGP.default_funcs(numinputs[end]))
    println(f,"numiterations: ",numiterations)
    #println(f,"nodearity: ",nodearity)
    #println(f,"active_only: ",active_only)
    println(f,"gl_repetitions: ",gl_repetitions)
    println(f,"perm_heuristic: ",perm_heuristic)
    #println(f,"max_steps",maxsteps)
    CSV.write( f, df, append=true, writeheader=true )
  end
  #println(df)
  df
end

function run_env_evolve!( rr::env_result_type; maxints_for_degen::Int64, base::Float64=2.0 )
  nodearity = 2   # built-in default
  p = Parameters( numinputs=rr.numinputs, numoutputs=rr.numoutputs, numinteriors=rr.numints, numlevelsback=rr.levelsback )
  funcs = default_funcs(rr.numinputs) 
  #print_parameters( p )
  gl = rand_env_goallist(rr.ngoals,rr.numinputs,rr.numoutputs,rr.gl_reps,rr.num_flip_bits,perturb_goal=rr.perturb_goal)
  #println("gl: ",gl)
  c = random_chromosome( p, funcs )
  #(rr.steps,rr.worse,rr.same,rr.better,c,output,goals_matched,matched_goals_list) = 
  #    mut_evolve(c,gl,funcs,rr.maxsteps,hamming_sel=rr.hamming_sel,active_only=rr.active_only)
  (c,rr.steps,rr.worse,rr.same,rr.better,output,matched_goals,matched_goals_list) = 
      mut_evolve(c,gl,funcs,rr.maxsteps,avgfitness=rr.avgfitness,perm_heuristic=rr.perm_heuristic,hamming_sel=rr.hamming_sel,fit_limit=rr.fit_limit)
  rr.avgfit = c.fitness
  rr.nactive = number_active( c )
  rr.redundancy = redundancy( c, base=base )
  rr.complexity = rr.numints <= maxints_for_degen ? complexity5( c, base=base ) : 0.0
  rr.degeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base ) : 0.0
  rr.sdegeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base, mutinf=mutinf2 ) : 0.0
  rr
end

function env_result( p::Parameters, num_goals::Int64, hamming_sel::Bool, active_only::Bool, max_steps::Int64, 
    gl_reps::Int64, num_flip_bits::Int64, perm_heuristic::Bool, perturb_goal::Bool, avgfitness::Bool, fit_limit::Float64 )
  env_result_type(
    p.numinputs,
    p.numoutputs,
    p.numinteriors,
    p.numlevelsback,
    num_goals,
    hamming_sel,
    active_only,
    max_steps,
    gl_reps,
    num_flip_bits,
    perturb_goal,
    avgfitness,
    perm_heuristic,
    fit_limit,
    0,      # steps
    0.0,    # avgfit
    0,      # same
    0,      # worse
    0,      # better
    0,      # nactive 
    0.0,    # redundancy
    0.0,    # complexity
    0.0,    # degeneracy
    0.0     # sdegeneract
  )
end

function env_result_to_tuple( rr::env_result_type )
 (
  rr.numinputs,
  rr.numoutputs,
  rr.numints,
  rr.levelsback,
  rr.ngoals,
  #rr.hamming_sel,
  #rr.active_only,
  rr.maxsteps,
  rr.gl_reps,
  rr.num_flip_bits,
  rr.perturb_goal,
  rr.avgfitness,
  rr.fit_limit,
  rr.steps,
  rr.avgfit,
  #rr.same,
  #rr.worse,
  #rr.better,
  rr.nactive,
  rr.redundancy,
  rr.complexity,
  rr.degeneracy,
  rr.sdegeneracy
 )
end
