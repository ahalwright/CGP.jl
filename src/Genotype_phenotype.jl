# Implement code for proposal evotech/complexity/cce/'Robustness evolvability compleity proposal 8_25_20.docx'.
# Includes evolving multiple chromosomes that correspond to every possible 3-input 1-output phenotype. 
# In other words, evolves 3-input 1-output chromosomes (genotypes) whose output is every possible 8-bit UInt.
# There are 
export run_geno_pheno_evolution, gp_result, gp_result_to_tuple, run_gp_evolve!
using DataFrames 
using CSV
using Statistics
using Distributions
geno_pheno_result_type = Main.CGP.geno_pheno_result_type
#context = construct_contexts(numinputs)[numinputs]
#p = Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
#iterations = 20
#print_parameters(p)
using CSV
using Dates

#function timer( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, csvfile )
function run_geno_pheno_evolution( iterations::Int64, numinputs::IntRange, numoutputs::IntRange, 
    numinteriors::IntRange, goallistlength::IntRange, maxsteps::IntRange,
    levelsback::IntRange, csvfile::String; 
    base::Float64=2.0, allgoals::Bool=false, active_only::Bool=false, gl_repetitions::IntRange=1)
  println(default_funcs(2))
  result =
    @timed run_geno_pheno( iterations, numinputs, numoutputs, numinteriors, goallistlength, maxsteps, levelsback,
        base=base, allgoals=allgoals, active_only=active_only, gl_repetitions=gl_repetitions )
  df = result[1]
  ttime = result[2]
  println("run time in minutes: ",ttime/60)
  println("cor( df.logsteps, df.complexity): ",cor( df.logsteps, df.complexity))
  println("cor( df.logsteps, df.gb_complexity): ",cor( df.logsteps, df.gb_complexity))
  println("cor( df.logsteps, df.degeneracy): ",cor( df.logsteps, df.degeneracy))
  println("cor( df.logsteps, df.sdegeneracy): ",cor( df.logsteps, df.sdegeneracy))
  println("cor( df.logsteps, df.f_mutinf): ",cor( df.logsteps, df.f_mutinf))
  println("cor( df.complexity, df.f_mutinf): ",cor( df.complexity, df.f_mutinf))
  println("cor( df.gb_complexity, df.f_mutinf): ",cor( df.gb_complexity, df.f_mutinf))
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  open( csvfile, "w" ) do f
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    println(f,"# run time in minutes: ",ttime/60)
    println(f,"# funcs: ", Main.CGP.default_funcs(numinputs[end]))
    println(f,"cor( df.logsteps, df.complexity): ",cor( df.logsteps, df.complexity))
    println(f,"cor( df.logsteps, df.gb_complexity): ",cor( df.logsteps, df.gb_complexity))
    println(f,"cor( df.logsteps, df.degeneracy): ",cor( df.logsteps, df.degeneracy))
    println(f,"cor( df.logsteps, df.sdegeneracy): ",cor( df.logsteps, df.sdegeneracy))
    println(f,"cor( df.logsteps, df.f_mutinf): ",cor( df.logsteps, df.f_mutinf))
    println(f,"cor( df.complexity, df.f_mutinf): ",cor( df.logsteps, df.f_mutinf))
    println(f,"cor( df.gb_complexity, df.f_mutinf): ",cor( df.gb_complexity, df.f_mutinf))
    #println(f,"# nodearity: ",nodearity)
    #println(f,"# active_only: ",active_only)
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

function run_geno_pheno( iterations::Int64, numinputs::IntRange, numoutputs::IntRange, 
    numinteriors::IntRange, goallistlength::IntRange, maxsteps::IntRange, levelsback::IntRange ; 
    base::Float64=2.0, allgoals::Bool=false, active_only::Bool=false, gl_repetitions::IntRange=1)
  max_numinteriors = collect(numinteriors)[end]
  test_MyInt(max_numinteriors)
  nodearity = 2
  hamming_sel = true
  all_gl = [ [[x]] for x = MyInt(0):MyInt(2^2^numinputs-1) ]
  repeat_limit=10
  gp_result_list = geno_pheno_result_type[]
  df = DataFrame() 
  df.goallist=Vector{Vector{MyInt}}[]
  df.numinputs=Int64[]
  df.numoutputs=Int64[]
  df.numints=Int64[]
  df.levelsback=Int64[]
  df.ngoals=Int64[]
  df.maxsteps=Int64[]
  df.gl_reps=Int64[]
  df.steps=Int64[]
  df.logsteps=Float64[]
  df.avgfit=Float64[]
  df.nactive=Int64[]
  df.redundancy=Float64[]
  df.complexity=Float64[]
  df.gb_complexity=Float64[]
  df.degeneracy=Float64[]
  df.sdegeneracy=Float64[]
  df.f_mutinf=Float64[]
  #println("size(df): ",size(df))
  for num_inputs = numinputs
    for num_outputs = numoutputs
      #fit_limit = Float64(num_outputs)
      funcs = default_funcs(num_inputs) 
      for num_goals = goallistlength
        for gl_reps = gl_repetitions 
          #gl = Vector{MyInt}[]
          for num_interiors = numinteriors
            for levsback = levelsback
              println("numinputs: ",num_inputs,"  numoutputs: ",num_outputs,"  numints: ",num_interiors,
                "  levsback: ", levsback, "  numgoals: ", num_goals)
              for max_steps = maxsteps
                for gl = (allgoals ? all_gl : [Vector{MyInt}[]])
                  for _ = 1:iterations
                    p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, levsback )
                    rr = gp_result( gl, p, num_goals, hamming_sel, active_only, max_steps, gl_reps )
                    push!(gp_result_list,rr)
                  end
                end
              end
            end
          end
        end
      end
    end
  end
  new_gp_result_list = pmap(r->run_gp_evolve!(r,maxints_for_degen=maxints_for_degen,gl_repetitions=gl_repetitions,base=base,repeat_limit=repeat_limit),gp_result_list)
  #new_gp_result_list = map(r->run_gp_evolve!(r,maxints_for_degen=Main.CGP.maxints_for_degen,gl_repetitions=gl_repetitions,base=base,repeat_limit=repeat_limit),gp_result_list)
  for r = new_gp_result_list
    new_row = gp_result_to_tuple(r)
    #println("length new row: ",length(new_row))
    Base.push!( df, new_row )
  end
  df
end

function run_gp_evolve!( rr::geno_pheno_result_type; maxints_for_degen::Int64, gl_repetitions::Int64=1, base::Float64=2.0,
      repeat_limit::Int64=10 )
  #repeat_limit = 10
  #println("run_gp_evolve! fault_tol: ",rr.fault_tol)
  nodearity = 2   # built-in default
  p = Parameters( numinputs=rr.numinputs, numoutputs=rr.numoutputs, numinteriors=rr.numints, numlevelsback=rr.levelsback )
  #print_parameters( p )
  if length(rr.goallist) > 0 
    println("rr.goallist: ",rr.goallist)
  end
  if length(rr.goallist) == 0
    rr.goallist = randgoallist(rr.ngoals,rr.numinputs,rr.numoutputs,repetitions=gl_repetitions)
  end
  funcs = default_funcs(rr.numinputs) 
  c = random_chromosome( p, funcs )
  repeat = 0
  (c,rr.steps,output,matched_goals,matched_goals_list) = mut_evolve(c,rr.goallist,funcs,rr.maxsteps)
  if rr.steps == rr.maxsteps
    println("goal: ",rr.goallist[1])
  end
  while repeat < repeat_limit && rr.steps == rr.maxsteps
    rr.numints += 1
    rr.levelsback += 1
    p_repeat = Parameters( numinputs=rr.numinputs, numoutputs=rr.numoutputs, numinteriors=rr.numints, numlevelsback=rr.levelsback )
    println("repeating function mut_evolve with numints: ",rr.numints,"  and with levsback: ",rr.levelsback)
    c = random_chromosome( p_repeat, funcs )
    (c,rr.steps,output,matched_goals,matched_goals_list) = 
        mut_evolve(c,rr.goallist,funcs,rr.maxsteps)
    #print_build_chromosome( c )
    repeat += 1
  end
  if repeat == repeat_limit
    println("run_gp_evolve!  repeat == repeat_limit ")
  end
  #if rr.numints > p.numinteriors
  #  println("rr.numints - p.numinteriors: ",rr.numints - p.numinteriors)
  #end
  rr.logsteps = log( 1.0 + rr.steps + rr.maxsteps*(rr.numints - p.numinteriors) )
  rr.avgfit = c.fitness
  rr.nactive = number_active( c )
  rr.redundancy = redundancy( c, base=base )
  rr.complexity = rr.numints <= maxints_for_degen ? complexity5( c, base=base ) : 0.0
  rr.gb_complexity = rr.numints <= maxints_for_degen ? gb_complexity_chrome( c, base=base ) : 0.0
  rr.degeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base ) : 0.0
  rr.sdegeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base, mutinf=mutinf2 ) : 0.0
  rr.f_mutinf = fmi_chrome( c )
  rr
end

function gp_result( gl::Vector{Vector{MyInt}}, p::Parameters, num_goals::Int64, hamming_sel::Bool, active_only::Bool, max_steps::Int64, 
      gl_reps::Int64 )
  geno_pheno_result_type(
    #Vector{MyInt}[],
    gl,
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
    0.0,    # logsteps
    0.0,    # avgfit
    0,      # nactive 
    0.0,    # redundancy
    0.0,    # complexity
    0.0,    # gb_complexity
    0.0,    # degeneracy
    0.0,    # sdegeneracy
    0.0     # f_mutinf
  )
end

function gp_result_to_tuple( rr::geno_pheno_result_type )
 (
  rr.goallist,
  rr.numinputs,
  rr.numoutputs,
  rr.numints,
  rr.levelsback,
  rr.ngoals,
  rr.maxsteps,
  rr.gl_reps,
  rr.steps,
  rr.logsteps,
  rr.avgfit,
  rr.nactive,
  rr.redundancy,
  rr.complexity,
  rr.gb_complexity,
  rr.degeneracy,
  rr.sdegeneracy,
  rr.f_mutinf
 )
end
