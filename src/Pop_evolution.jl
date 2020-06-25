# Implements population based evoluiont with runs of individual evolution in each population generation
# Populations are lists of chromosomes
export pop_evolve, run_pop_evolve_df, run_pop_evolve!, pop_result, propsel!, tournsel!
using DataFrames 
using CSV
using Statistics
using Distributions
using Distributed
pop_result_type = Main.CGP.pop_result_type
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
max_pop_gens_rng = 10:10
max_indiv_steps_rng = 20:20
hamming_rng = true:true
robust_rng = true:true
run_pop_evolve_df( iterations, numinputs, numoutputs, numinteriors, ngoals, levelsback, 
    max_pop_gens_rng, max_indiv_steps_rng, hamming_rng, robust_rng, "testdata.csv" )
=#
#context = construct_contexts(numinputs)[numinputs]
#p = Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
#iterations = 20
#print_parameters(p)

function run_pop_evolve_df( numiterations::Int64, numinputs::AbstractRange{Int64}, numoutputs::AbstractRange{Int64}, 
    numinteriors::AbstractRange{Int64}, ngoals_rng::AbstractRange{Int64}, 
    levelsback::AbstractRange, max_pop_gens_rng::AbstractRange{Int64}, max_indiv_steps_rng::AbstractRange{Int64},
    popsize_rng::AbstractRange{Int64},
    hamming_rng::AbstractRange{Bool}, robust_rng::AbstractRange{Bool}, csvfile::String; 
    base::Float64=2.0, active_only::Bool=false, gl_repetitions::Int64=1 )
  maxints_for_degen = 20
  all_max_sel=true
  nodearity = 2
  pop_result_list = pop_result_type[]
  df = DataFrame() 
  df.numinputs=Int64[]
  df.numoutputs=Int64[]
  df.numints=Int64[]
  df.levelsback=Int64[]
  df.ngoals=Int64[]
  df.max_indiv_steps=Int64[]
  df.max_pop_gens=Int64[]
  df.popsize=Int64[]
  df.robust_sel=Bool[]
  df.all_max_sel = Bool[]
  #df.active_only=Bool[]
  df.maxfit=Float64[]
  df.maxrobust=Float64[]
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
        for num_goals = ngoals_rng
          println("numinputs: ",num_inputs,"  numoutputs: ",num_outputs,"  numints: ",num_interiors,"  numgoals: ",num_goals)
          for levsback = levelsback
            for max_pop_gens = max_pop_gens_rng
              for max_indiv_steps = max_indiv_steps_rng
                for popsize = popsize_rng
                  for hamming_sel = hamming_rng
                    for robust_sel = robust_rng
                      for active_only = [false]
                        println("hamming_sel: ",hamming_sel,"  robust_sel: ",robust_sel)
                        for _ = 1:numiterations
                          p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, levsback )
                          rr = pop_result( p, ngoals=num_goals, max_pop_gens=max_pop_gens, popsize=popsize,
                              max_indiv_steps=max_indiv_steps, hamming_sel=hamming_sel, robust_sel=robust_sel, 
                              all_max_sel=all_max_sel, active_only=active_only )
                          push!(pop_result_list,rr)
                          #new_row = pop_result_to_tuple(rr)
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
  new_pop_result_list = pmap(r->run_pop_evolve!(r,maxints_for_degen=maxints_for_degen,gl_repetitions=gl_repetitions,base=base),pop_result_list)
  #new_pop_result_list = map(r->run_pop_evolve!(r,maxints_for_degen=maxints_for_degen,gl_repetitions=gl_repetitions,base=base),pop_result_list)
  for r = new_pop_result_list
    new_row = pop_result_to_tuple(r)
    Base.push!( df, new_row )
  end
  println(default_funcs(2))
  open( csvfile, "w" ) do f
    println(f,"MyInt: ",Main.CGP.MyInt)
    println(f,"funcs: ", Main.CGP.default_funcs(numinputs[end]))
    println(f,"nodearity: ",nodearity)
    println(f,"gl_repetitions: ",gl_repetitions)
    println(f,"maxints_for_degen: ",maxints_for_degen)
    println(f,"all_max_sel: ",all_max_sel)
    CSV.write( f, df, append=true, writeheader=true )
  end
  #println(df)
  df
end

function run_pop_evolve!( rr::pop_result_type; maxints_for_degen::Int64, gl_repetitions::Int64=1, base::Float64=2.0 )
  nodearity = 2   # built-in default
  p = Parameters( numinputs=rr.numinputs, numoutputs=rr.numoutputs, numinteriors=rr.numints, numlevelsback=rr.levelsback )
  #print_parameters( p )
  gl = randgoallist(rr.ngoals,rr.numinputs,rr.numoutputs,repetitions=gl_repetitions)
  #println("gl: ",gl)
  funcs = default_funcs(rr.numinputs) 
  c = random_chromosome( p, funcs )
  (rr.maxfit, maxfit_c, rr.maxrobust, maxrobust_c, pop ) =
      pop_evolve(p,rr.popsize,rr.max_pop_gens,rr.max_indiv_steps,gl,funcs,
      hamming_sel=rr.hamming_sel,robust_sel=rr.robust_sel,all_max_sel=rr.all_max_sel)
  if rr.robust_sel
    c = maxrobust_c
  else
    c = maxfit_c
  end
  rr.nactive = number_active( c )
  rr.redundancy = redundancy( c, base=base )
  rr.complexity = rr.numints <= maxints_for_degen ? complexity5( c, base=base ) : 0.0
  rr.degeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base ) : 0.0
  rr.sdegeneracy = rr.numints <= maxints_for_degen ? degeneracy( c, base=base, mutinf=mutinf2 ) : 0.0
  rr
end

function pop_evolve( params::Parameters, popsize::Int64, max_pop_gens::Int64, max_indiv_steps::Int64,
    goallist::GoalList, funcs::Vector{Func}; uniform_start::Bool=false, hamming_sel::Bool=true,
    robust_sel::Bool=true, all_max_sel::Bool=false, active_only::Bool=false )
  println("pop evolve: hamming_sel: ",hamming_sel, "  robust_sel: ",robust_sel,"  all_max_sel: ",all_max_sel)
  if uniform_start
    c = random_chromosome(params,funcs)
    pop = fill( c, popsize )
  else
    pop = fill( random_chromosome(params,funcs), popsize )
  end
  max_gen_fitness = zeros(Float64,max_pop_gens)
  max_gen_fit_chrome = Chromosome[]
  max_gen_robustness = zeros(Float64,max_pop_gens)
  max_gen_robust_chrome = Chromosome[]
  for g = 1:max_pop_gens
    using_robustness = false  # Note if all_max_sel is not true, the maximum fitness individual may be lost in robustness selection
    println("gen: ",g,"  fitnesses: ",[pop[i].fitness for i =1:popsize]," robustnesses: ",[Base.round(pop[i].robustness,digits=2) for i =1:popsize])
    for i = 1:popsize
      (pop[i],step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( pop[i], goallist, funcs, max_indiv_steps, hamming_sel=hamming_sel )
      println("(g,i): ",(g,i),"  pop[i].fitness: ",pop[i].fitness,"  pop[i].robustness: ",pop[i].robustness)
    end
    # Do proportional selection
    fitness_vector = [ pop[i].fitness for i = 1:popsize ]
    (maxfit, maxfit_indices ) = findmaxall(fitness_vector)
    #println("gen: ",g,"  max fit before selection: ",maxfit,"  pop indices: ",maxfit_indices)
    if using_robustness
      robustness_vector = [ pop[i].robustness for i = 1:popsize ]
      (maxrobust, maxrobust_indices ) = findmaxall(robustness_vector)
      #println("gen: ",g,"  max robust before selection: ",maxrobust,"  pop indices: ",maxrobust_indices)
    end
    if (!robust_sel && !all_max_sel && maxfit == params.numoutputs)
      println("found single individual optimum so break out of generation loop")
      fitness_vector = [ pop[i].fitness for i = 1:popsize ]
      (max_gen_fitness[g],maxfit_indices) = findmaxall(fitness_vector)
      push!(max_gen_fit_chrome,pop[maxfit_indices[1]])
      break  # found single individual optimum so break out of loop
    elseif (!robust_sel && all_max_sel && maxfit == params.numoutputs && length(maxfit_indices) == popsize)
      println("found optimum in entire population so break out of generation loop")
      fitness_vector = [ pop[i].fitness for i = 1:popsize ]
      (max_gen_fitness[g],maxfit_indices) = findmaxall(fitness_vector)
      push!(max_gen_fit_chrome,pop[maxfit_indices[1]])
      break  # found optimum in entire population so break out of loop
    elseif robust_sel && maxfit == params.numoutputs && ( !all_max_sel || (all_max_sel && length(maxfit_indices) == popsize))
      println("doing robustness selection")
      using_robustness = true
      for i = 1:popsize
        pop[i].robustness = mutational_robustness( pop[i], funcs, active_only=active_only )
      end
      robustness_vector = [ pop[i].robustness for i = 1:popsize ]
      propsel!( pop, robustness_vector, maxfit=findmax(robustness_vector)[1] ) 
    else
      println("propsel")
      propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] ) 
    end 
    fitness_vector = [ pop[i].fitness for i = 1:popsize ]
    (max_gen_fitness[g],maxfit_indices) = findmaxall(fitness_vector)
    push!(max_gen_fit_chrome,pop[maxfit_indices[1]])
    println("gen: ",g,"  max fit after selection: ",max_gen_fitness[g],"  at index: ",maxfit_indices)
    robustness_vector = [ pop[i].robustness for i = 1:popsize ]
    (max_gen_robustness[g],maxrobust_indices) = findmaxall(robustness_vector)
    push!(max_gen_robust_chrome,pop[maxrobust_indices[1]])
    #push!(max_gen_fit_chrome, pop[maxfit_indices[1]])
    if using_robustness
      println("gen: ",g,"  max robust after selection: ",max_gen_robustness[g],"  at index: ",maxrobust_indices)
    end
  end
  println("length(max_gen_fit_chrome): ",length(max_gen_fit_chrome))
  println("length(max_gen_robust_chrome): ",length(max_gen_robust_chrome))
  (max_fit_all_gens, max_fit_all_gens_ind) = findmax(max_gen_fitness)
  #push!(max_gen_fit_chrome, max_gen_fit_chrome[max_fit_all_gens_ind])
  (max_robust_all_gens, max_robust_all_gens_ind) = findmax(max_gen_robustness)
  println("max fit all gens: ",max_fit_all_gens,"  max robust all gens: ",max_robust_all_gens)
  maxfit_gen_chromosome = length(max_gen_fit_chrome) > 0 ? max_gen_fit_chrome[max_fit_all_gens_ind] : nothing
  maxrobust_gen_chromosome = length(max_gen_robust_chrome) > 0 ? max_gen_robust_chrome[max_robust_all_gens_ind] : nothing
  return (max_fit_all_gens,maxfit_gen_chromosome,
    max_robust_all_gens,maxrobust_gen_chromosome ,pop)
end

function pop_result( ; numinputs::Int64=2, numoutputs::Int64=2, numints::Int64=6, levelsback::Int64=6,
    ngoals::Int64=4, popsize::Int64=5, max_pop_gens::Int64=10, max_indiv_steps::Int64=20, 
    hamming_sel::Bool=true, robust_sel::Bool=true, all_max_sel::Bool=true, active_only::Bool=false,
    maxfit::Float64=0.0, maxrobust::Float64=0.0, nactive::Int64=0, redundancy::Float64=0.0, 
    complexity::Float64=0.0, degeneracy::Float64=0.0, sdegeneracy::Float64=0.0)
  pop_result_type( numinputs, numoutputs, numints, levelsback, ngoals, popsize, max_pop_gens, max_indiv_steps,
      hamming_sel, robust_sel, all_max_sel, active_only, maxfit, maxrobust, nactive, redundancy, complexity,
      degeneracy, sdegeneracy )
end

function pop_result( p::Parameters ; 
    ngoals::Int64=4, popsize::Int64=3, max_pop_gens::Int64=10, max_indiv_steps::Int64=20, 
    hamming_sel::Bool=true, robust_sel::Bool=true, all_max_sel::Bool=true, active_only::Bool=false,
    maxfit::Float64=0.0, maxrobust::Float64=0.0, nactive::Int64=0, redundancy::Float64=0.0, complexity::Float64=0.0, 
    degeneracy::Float64=0.0, sdegeneracy::Float64=0.0)
    pop_result_type( p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, ngoals, popsize, max_pop_gens, 
    max_indiv_steps, hamming_sel, robust_sel, all_max_sel, active_only, maxfit, maxrobust, nactive, redundancy, 
    complexity, degeneracy, sdegeneracy )
end

function pop_result_to_tuple( rr::pop_result_type )
  # must match with the columns of the dataframe defined at the beginning of function run_pop_evolve_df()
  (
    rr.numinputs,
    rr.numoutputs,
    rr.numints,
    rr.levelsback,
    rr.ngoals,
    rr.max_indiv_steps,
    rr.max_pop_gens,
    rr.popsize,
    rr.robust_sel,
    rr.all_max_sel,
    #rr.active_only,
    rr.maxfit,
    rr.maxrobust,
    rr.nactive,
    rr.redundancy,
    rr.complexity,
    rr.degeneracy,
    rr.sdegeneracy,
  )
end

@doc """function propsel!(p::Population, fitness::Vector{Float64} )
Conduct proportional selection in-place.
"""
function propsel!( pop::Vector{Chromosome}, fitness::Vector{Float64}; maxfit::Float64=0.0  )
  if maxfit == 0.0
    maxfit = maximum(fitness)
    if maxfit == 0.0
      # all elements have fitness zero
      return
    end
  end
  n = length(pop)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = pop[i].fitness / maxfit
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end 

function tournsel!( pop::Vector{Chromosome}, fitness::Vector{Float64}, tournsize::Int64=2 )
  popsize = length(pop)
  newpop = Chromosome[]
  for k = 1:popsize
    indices = rand(1:popsize,tournsize)
    #println("indices: ",indices)
    (max,index) = findmax( [(pop[i].fitness,i) for i in indices] )
    #println("(max,index): ",(max,index))
    push!(newpop,deepcopy(pop[max[2]]))
  end
  newpop
end
