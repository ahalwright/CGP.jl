# An evolution that first does mutational evolution using goal fitness and then 
#    population evolultion using mutational robustness as fitness

export mut_pop_evolve, run_mut_pop_evolve!, run_mut_pop_evolve_df, round2, round3

using DataFrames
using CSV
using Statistics
using Distributions
using Distributed

# Note that selection on robustness is always done, and population evolution is not started until
#   popsize individuals have been found that achieve maximum goal fitness.
# Thus, pop_gens should be set to a very large value
function mut_pop_evolve( nreps::Int64, params::Parameters, popsize::Int64, pop_gens::Int64, indiv_steps::Int64,
    goallist::GoalList, funcs::Vector{Func}, csvfile::String; tourn_size::Int64=0,  # 0 means use prosel
    uniform_start::Bool=false, hamming_sel::Bool=true, active_only::Bool=false )
  println("Int(trunc(indiv_steps/10)): ",Int(trunc(indiv_steps/10)))
  popsize_multiplier = 2
  initial_popsize = popsize_multiplier*popsize 
  found_optimum = fill(false,initial_popsize)
  null_output = fill(MyInt(0),params.numoutputs)
  output_list = fill(null_output,initial_popsize)
  if uniform_start
    c = random_chromosome(params,funcs)
    pop = fill( c, initial_popsize )
  else
    pop = [ random_chromosome(params,funcs) for _=1:initial_popsize]
  end
  count_at_optimal_fitness = 0
  mut_step = 0
  max_mut_steps = 500  # prevent an infinite loop
  while mut_step < max_mut_steps && count_at_optimal_fitness < popsize
    for i = 1:initial_popsize
      if !found_optimum[i]
        (pop[i],step,worse,same,better,output_list[i],matched_goals,matched_goals_list) = 
            mut_evolve( pop[i], goallist, funcs, Int(trunc(indiv_steps/10)))
        found_optimum[i] = (sort(output_list[i]) == sort(goallist[matched_goals[1]]))
        count_at_optimal_fitness += found_optimum[i]   # in Julia, true==1 and false==0
      end
    end
    mut_step += 1
    println("mut_step: ",mut_step,"  count_at_optimal_fitness: ",count_at_optimal_fitness)
  end    
  if mut_step == max_mut_steps
    error(" mut_step == max_mut_steps in mut_pop_evolve" )
  end
  #println("found_optimum: ",found_optimum)
  npop = pop[ found_optimum ][1:popsize]
  #println("count_at_optimal: ",count_at_optimal_fitness,"  mut_step: ",mut_step )
  #(count_at_optimal_fitness, mut_step, npop  )
  #mean_robustness = zeros(Float64,pop_gens)
  #max_robustness = zeros(Float64,pop_gens)
  df = DataFrame()
  df.generation = collect(1:pop_gens)
  df.numinputs = fill(params.numinputs,pop_gens)
  df.numoutputs = fill(params.numoutputs,pop_gens)
  df.numints = fill(params.numinteriors,pop_gens)
  df.levsback = fill(params.numlevelsback,pop_gens)
  df.mean_robustness = zeros(Float64,pop_gens)
  df.max_robustness = zeros(Float64,pop_gens)
  df.nactive = zeros(Float64,pop_gens)
  df.complexity = zeros(Float64,pop_gens)
  df.degeneracy = zeros(Float64,pop_gens)
  for r = 1:nreps
    println("nrep: ",r)
    for g = 1:pop_gens
      #println("gen: ",g)
      for i = 1:popsize
        prev_output = output_values(npop[i])
        if g > 1
          #print("bef na: ",number_active(npop[i]))
          npop[i] = neutral_mutate( npop[i], funcs ) 
          #println("  aft na: ",number_active(npop[i]))
          #mutate_chromosome!(npop[i],funcs)
        end
        npop[i].robustness = mutational_robustness( npop[i], funcs, active_only=active_only )
        new_output = output_values(npop[i])
        if prev_output != new_output
          println("prevout: ",prev_output,"  newout:  ",new_output)
        end
      end
      robustness_vector = [ npop[i].robustness for i = 1:popsize ]
      nactive_vector = [ number_active( npop[i] ) for i = 1:popsize ]
      complexity_vector = [ complexity5( npop[i] ) for i = 1:popsize ]
      degeneracy_vector = [ degeneracy( npop[i] ) for i = 1:popsize ]
      df.mean_robustness[g] += mean(robustness_vector)
      df.max_robustness[g] += maximum(robustness_vector)
      df.nactive[g] += mean(nactive_vector)
      df.complexity[g] += mean(complexity_vector)
      df.degeneracy[g] += mean(degeneracy_vector)
      #print("  rv: ",[round2(rv) for rv in robustness_vector],"  ")
      #println("  na: ",[na for na in nactive_vector],"  ")
      #println( "mean: ",round2(mean_robustness[g]),"  max: ",round2(max_robustness[g]))
      if tourn_size == 0
        propsel!( npop, robustness_vector, maxfit=df.max_robustness[g])
      else
        tournsel!( npop, robustness_vector, tourn_size)
      end
    end
  end
  df.mean_robustness /= nreps
  df.max_robustness /= nreps
  df.nactive /= nreps
  df.complexity /= nreps
  df.degeneracy /= nreps
  #=
  for g = 1:pop_gens
    print("g:",g,"  ")
    println( "mean: ",round3(df.mean_robustness[g]),"  max: ",round3(df.max_robustness[g]),
        "  na: ", round3(df.nactive[g]))
  end
  =#
  open( csvfile, "w" ) do f
    println(f,"MyInt: ",Main.CGP.MyInt)
    println(f,"funcs: ", Main.CGP.default_funcs(params.numinputs))
    println(f,"nodearity: ",params.nodearity)
    println(f,"nreps: ",nreps)
    println(f,"popsize: ",popsize)
    println(f,"tourn_size: ", tourn_size)
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

# Mutates chromosome c until a chromosome with the same output is created.
# The returned chromosome may be multiple mutations away from c.
function neutral_mutate( c::Chromosome, funcs::Vector{Func}  )
  orig_c = deepcopy(c)
  orig_output = output_values(c)
  #println("orig output: ",orig_output)
  mutate_chromosome!( c, funcs )
  output = output_values(c)
  #println("bef output: ",output)
  count = 0
  while count < 200 &&  output != orig_output
    c = deepcopy(orig_c)
    mutate_chromosome!( c, funcs )
    output = output_values(c)
    #println("lll output: ",output)
    count += 1
  end
  c
end
    

round2( x::Float64 ) = Base.round(x,digits=2)
round3( x::Float64 ) = Base.round(x,digits=3)
