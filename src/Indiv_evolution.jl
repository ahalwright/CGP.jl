# Implement experiments to measure the success of subgoal evolution as a function of:
#  1.  Funcs
#  2.  numinteriors
#  3.  numinputs #  4.  numoutputs
#  5.  Length of goallist
#  5.  max_steps
# Uses randomly generated goallist
# Keep track of number of "worse" and "same" updates of Chromosome since its last improvement
export run_mut_evolution
using DataFrames 
using CSV
using Statistics
using Distributions
indiv_result_type = Main.CGP.indiv_result_type
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

function run_mut_evolution( numiterations::Int64, numinputs::IntRange, numoutputs::IntRange, 
    numinteriors::IntRange, goallistlength::IntRange, maxsteps::IntRange,
    levelsback::IntRange, hamming_rng::IntRange, csvfile::String; 
    base::Float64=2.0, active_only::Bool=false, gl_repetitions::IntRange=1 )
  maxints_for_degen = 20
  nodearity = 2
  run_result_list = indiv_result_type[]
  df = DataFrame() 
  df.numinputs=Int64[]
  df.numoutputs=Int64[]
  df.numints=Int64[]
  df.levelsback=Int64[]
  df.ngoals=Int64[]
  df.hamming_sel=Bool[]
  #df.robust_sel=Bool[]
  df.active_only=Bool[]
  df.maxsteps=Int64[]
  df.gl_reps=Int64[]
  df.steps=Int64[]
  df.same=Int64[]
  df.worse=Int64[]
  df.better=Int64[]
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
            for hamming_sel = hamming_rng
                for active_only = [false]
                  for gl_reps = gl_repetitions 
                    println("hamming_sel: ",hamming_sel,"  active_only: ",active_only)
                    for max_steps = maxsteps
                      for _ = 1:numiterations
                        p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, levsback )
                        rr = run_result( p, num_goals, hamming_sel, active_only, max_steps, gl_reps )
                        push!(run_result_list,rr)
                        #run_mut_evolve!(rr)
                        #new_row = run_result_to_tuple(rr)
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
  new_run_result_list = pmap(r->run_mut_evolve!(r,maxints_for_degen=maxints_for_degen,gl_repetitions=gl_repetitions,base=base),run_result_list)
  #new_run_result_list = map(r->run_mut_evolve!(r,maxints_for_degen=maxints_for_degen,gl_repetitions=gl_repetitions,base=base),run_result_list)
  for r = new_run_result_list
     new_row = run_result_to_tuple(r)
     Base.push!( df, new_row )
  end
  println(default_funcs(2))
  open( csvfile, "w" ) do f
    println(f,"funcs: ", Main.CGP.default_funcs(numinputs[end]))
    println(f,"nodearity: ",nodearity)
    #println(f,"active_only: ",active_only)
    println(f,"gl_repetitions: ",gl_repetitions)
    println(f,"max_steps",maxsteps)
    CSV.write( f, df, append=true, writeheader=true )
  end
  #println(df)
  df
end

function run_mut_evolve!( rr::indiv_result_type; maxints_for_degen::Int64, gl_repetitions::Int64=1, base::Float64=2.0 )
  nodearity = 2   # built-in default
  p = Parameters( numinputs=rr.numinputs, numoutputs=rr.numoutputs, numinteriors=rr.numints, numlevelsback=rr.levelsback )
  #print_parameters( p )
  gl = randgoallist(rr.ngoals,rr.numinputs,rr.numoutputs,repetitions=gl_repetitions)
  #println("gl: ",gl)
  funcs = default_funcs(rr.numinputs) 
  c = random_chromosome( p, funcs )
  #(rr.steps,rr.worse,rr.same,rr.better,c,output,goals_matched,matched_goals_list) = mut_evolve(c,gl,funcs,rr.maxsteps,hamming_sel=rr.hamming_sel,active_only=rr.active_only)
  (c,rr.steps,rr.worse,rr.same,rr.better,output,matched_goals,matched_goals_list) = mut_evolve(c,gl,funcs,rr.maxsteps,hamming_sel=rr.hamming_sel)
  rr.nactive = number_active( c )
  rr.redundancy = redundancy( c, base=base )
  rr.complexity = rr.numints <= maxints_for_degen ? complexity5( c, base=base ) : 0.0
  rr.degeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base ) : 0.0
  rr.sdegeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base, mutinf=mutinf2 ) : 0.0
  rr
end

function run_result( p::Parameters, num_goals::Int64, hamming_sel::Bool, active_only::Bool, max_steps::Int64, gl_reps::Int64 )
  indiv_result_type(
    p.numinputs,
    p.numoutputs,
    p.numinteriors,
    p.numlevelsback,
    num_goals,
    hamming_sel,
    active_only,
    max_steps,
    gl_reps,
    0,      # steps
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

function run_result_to_tuple( rr::indiv_result_type )
 (
  rr.numinputs,
  rr.numoutputs,
  rr.numints,
  rr.levelsback,
  rr.ngoals,
  rr.hamming_sel,
  rr.active_only,
  rr.maxsteps,
  rr.gl_reps,
  rr.steps,
  rr.same,
  rr.worse,
  rr.better,
  rr.nactive,
  rr.redundancy,
  rr.complexity,
  rr.degeneracy,
  rr.sdegeneracy
 )
end