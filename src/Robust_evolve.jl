# A version of evolution with selection for robustness that evolves population until all individuals are at optimal goal fitness,
#    then evolves with selection for mutational robustness for a fixed number of generations.
export run_robust_evolve, robust_evolve
using DataFrames
using Statistics
using CSV
using Distributed

# To run julia:  julia -i -p 4 -L CGP.jl

function run_robust_evolve( nreps::Int64, p::Parameters, popsize::Int64, max_pop_gens::Int64, indiv_steps::Int64, robust_steps::Int64, fit_steps::Int64,
      goallist::GoalList, tourn_size::Int64, csvfile::String )   # tourn_size==0 means proportional selection
  funcs = default_funcs(p.numinputs)
  test_MyInt(p.numinteriors) 
  df_results = DataFrame[]
  df = DataFrame()
  df.generation = collect(1:2*robust_steps)
  df.numinputs = fill(p.numinputs,2*robust_steps)
  df.numoutputs = fill(p.numoutputs,2*robust_steps)
  df.numints = fill(p.numinteriors,2*robust_steps)
  df.levsback = fill(p.numlevelsback,2*robust_steps)
  df.popsize = fill(popsize,2*robust_steps)
  df.robust_sel = fill(false,2*robust_steps)
  df.mean_robfit = zeros(Float64,2*robust_steps)
  df.max_robfit = zeros(Float64,2*robust_steps)
  df.mean_nactive = zeros(Float64,2*robust_steps)
  df.var_nactive = zeros(Float64,2*robust_steps)
  df.mean_complexity = zeros(Float64,2*robust_steps)
  df.var_complexity = zeros(Float64,2*robust_steps)
  df.mean_degeneracy = zeros(Float64,2*robust_steps)          
  df.var_degeneracy = zeros(Float64,2*robust_steps)          
  df_init = [ deepcopy(df) for _ = 1:nreps]
  df_results = pmap( ddf->robust_evolve(ddf,p, popsize, max_pop_gens, indiv_steps, robust_steps, fit_steps, goallist, tourn_size ), df_init)
  #df_results = map( ddf->robust_evolve(ddf,p, popsize, max_pop_gens, indiv_steps, robust_steps, fit_steps, goallist, tourn_size ), df_init)
  #println("df_results.mean_complexity: ",[df_results[i].mean_complexity for i = 1:nreps])
  #println("df_results.mean_complexity: ",[d.mean_complexity for d in df_results])
  df.robust_sel = df_results[1].robust_sel
  # Note that df.mean_robfit is a vector over reps of vectors over generations
  df.mean_robfit = mean([d.mean_robfit for d in df_results]) 
  df.max_robfit = mean([d.max_robfit for d in df_results])
  # Note that df.mean_nactive is a vector over reps of vectors over generations
  # the mean() function computes the mean over reps for each generation
  # the var() function computes the sample variance over reps for each generation
  df.mean_nactive = mean([d.mean_nactive for d in df_results])
  @assert df.mean_nactive == sum([d.mean_nactive for d in df_results])/nreps
  df.var_nactive = var([d.mean_nactive for d in df_results])
  df.mean_complexity = mean([d.mean_complexity for d in df_results])
  df.var_complexity = var([d.mean_complexity for d in df_results])
  #println("df.var_complexity: ",df.var_complexity)
  df.mean_degeneracy = mean([d.mean_degeneracy for d in df_results])
  df.var_degeneracy = var([d.mean_degeneracy for d in df_results])
  open( csvfile, "w" ) do f
    println(f,"MyInt: ",Main.CGP.MyInt)
    println(f,"funcs: ",funcs)
    println(f,"nreps: ",nreps)
    println(f,"indiv_steps: ",indiv_steps)
    println(f,"robust_steps: ",robust_steps)
    println(f,"fit_steps: ",fit_steps)
    println(f,"popsize: ",popsize)
    println(f,"max_pop_gens: ",max_pop_gens)
    println(f,"ngoals: ",length(goallist))
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

function robust_evolve( df::DataFrame, p::Parameters, popsize::Int64, max_pop_gens::Int64, indiv_steps::Int64, robust_steps::Int64, fit_steps::Int64,
    goallist::GoalList, tourn_size::Int64 )   # tourn_size==0 means proportional selection
  funcs = default_funcs(p.numinputs)
  for robust_sel = false:true
    pop = [ random_chromosome(p,funcs) for _ = 1:popsize ]
    gg = 1
    while gg < max_pop_gens
      #println("gen: ",gg,"  fitnesses: ",[pop[i].fitness for i =1:popsize] )
      for i = 1:popsize
        (pop[i],step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( pop[i], goallist, funcs, indiv_steps )
        #println("(gg,i): ",(gg,i),"  pop[i].fitness: ",pop[i].fitness)
      end
      fitness_vector = [ pop[i].fitness for i = 1:popsize ]
      if minimum(fitness_vector) == Float64(p.numoutputs)
        break
      end
      gg += 1
    end
    if gg == max_pop_gens
      error("generation_limit_exceeded")
    else
      println("population optimized at genertion ",gg)
      #println("fitvec: ", [ pop[i].fitness for i = 1:popsize ])
    end
    gr_increment = robust_sel*robust_steps
    #println("gr_increment: ",gr_increment)
    for gr = (gr_increment + 1):(gr_increment + robust_steps)
      df.robust_sel[gr] = robust_sel
      for i = 1:popsize
        if robust_sel
          pop[i].robustness = mutational_robustness( pop[i], funcs )
          robfit_vector = [ pop[i].robustness for i = 1:popsize ] 
        else
          (pop[i],step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( pop[i], goallist, funcs, fit_steps, hamming_sel=robust_sel) 
          robfit_vector = [ pop[i].fitness for i = 1:popsize ] 
        end
      end        
      #println("robfit before selection: ", [ round2(v) for v in robfit_vector] )
      if tourn_size == 0
        propsel!( pop, robfit_vector, maxfit=findmax(robfit_vector)[1] )
      else
        tournsel!( pop, robfit_vector, tourn_size )
      end
      robfit_vector = robust_sel ? [ pop[i].robustness for i = 1:popsize ] : [ pop[i].fitness for i = 1:popsize ] 
      nactive_vector = [ number_active( pop[i] ) for i = 1:popsize ]
      complexity_vector = df.numints <= maxints_for_degen ? [ complexity5( pop[i] ) for i = 1:popsize ] : zeros(MyInt,popsiz)
      degeneracy_vector = df.numints <= maxints_for_degen ? [ degeneracy( pop[i] ) for i = 1:popsize ] : zeros(MyInt,popsize)
      df.mean_robfit[gr] = mean(robfit_vector)
      df.max_robfit[gr] = maximum(robfit_vector)
      df.mean_nactive[gr] = mean(nactive_vector)
      df.mean_complexity[gr] = mean(complexity_vector)
      df.mean_degeneracy[gr] = mean(degeneracy_vector)    
      #println("robfit after selection: ", [ round2(v) for v in robfit_vector] )
      #=
      if gr==7 || gr==14
        println("gr: ",gr,"  complexity: ",complexity_vector,"  mean compleixty: ",mean(complexity_vector))
      end
      =#
    end
    #robfit_vector = [ pop[i].robfit for i = 1:popsize ]
  end
  #println("robust_evolve finished")
  df
end
