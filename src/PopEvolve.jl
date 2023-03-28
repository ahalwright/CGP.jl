# Objective:  Implement a population-based evolution for G-P maps.
using Distributions
# A population is a Vector of Chromosomes or LinCircuits.
  
function run_pop_evolve( p::Parameters, funcs::Vector{Func}, target::Goal, popsize::Int64, mutrate::Float64, ngens::Int64, nreps::Int64; 
      use_cgp::Bool=true, pop::Union{Vector{Chromosome},Vector{LinCircuit}}=Chromosome[] )
  meanfit_list = Float64[] 
  maxfit_list = Float64[] 
  first_opt_gen_list = Int64[]
  first_opt_majority_list = Int64[]
  for rep = 1:nreps
    ( meanfit, maxfit, first_opt_gen, first_majority_gen ) = pop_evolve( p, funcs, target, popsize, mutrate, ngens, use_cgp=use_cgp, pop=pop )
    push!(meanfit_list, meanfit )
    push!(maxfit_list, maxfit )
    push!(first_opt_gen_list, first_opt_gen )
    push!(first_opt_majority_list, first_majority_gen )
  end
  df = DataFrame( :target=>fill(target,nreps), 
                  :popsize=>fill(popsize,nreps), 
                  :mutrate=>fill(mutrate,nreps), 
                  :ngens=>fill(ngens,nreps), 
                  :meanfit=>meanfit_list,
                  :maxfit=>maxfit_list,
                  :first_opt_gen=>first_opt_gen_list,
                  :first_opt_majority=>first_opt_majority_list )
end

# Do nreps runs of pop_evolve() and return the fraction of runs where maxfit==1 and the fraction of runs where first_opt_majority>0
# nfunct_evals is the product of popsize and ngens.
function summary_pop_evolve( p::Parameters, funcs::Vector{Func}, target::Goal, mutrate_list::Vector{Float64}, nfunct_evals::Int64, popsize_list::Vector{Int64}, nreps::Int64; 
      use_cgp::Bool=true, pop::Union{Vector{Chromosome},Vector{LinCircuit}}=Chromosome[] )
  result_list = Tuple{Float64,Int64,Float64,Float64}[]
  for mutrate in mutrate_list
    for popsize  in popsize_list
      ngens = Int(ceil( nfunct_evals/popsize ))
      #println("mutrate: ",mutrate,"  popsize: ",popsize,"  ngens: ",ngens)
      meanfit_list = Float64[]
      maxfit_list = Float64[]
      first_opt_gen_list = Int64[]
      first_opt_majority_list = Int64[]
      for rep = 1:nreps
        ( meanfit, maxfit, first_opt_gen, first_majority_gen ) = pop_evolve( p, funcs, target, popsize, mutrate, ngens, use_cgp=use_cgp, pop=pop )
        #println("rep: ",rep,"  ( meanfit, maxfit, first_opt_gen, first_majority_gen ): ", ( meanfit, maxfit, first_opt_gen, first_majority_gen ) )
        push!(meanfit_list, meanfit )
        push!(maxfit_list, maxfit )
        push!(first_opt_gen_list, first_opt_gen )
        push!(first_opt_majority_list, first_majority_gen )
      end
      fract_maxfit = length(findall(x->x==1.0,maxfit_list)) / nreps
      fract_first_majority = length(findall(x->x>0,first_opt_majority_list)) / nreps
      push!( result_list, ( mutrate, popsize, fract_maxfit, fract_first_majority) )
    end
  end
  result_list
end

function pop_evolve( p::Parameters, funcs::Vector{Func}, target::Goal, popsize::Int64, mutrate::Float64, ngens::Int64; 
      use_cgp::Bool=true, pop::Union{Vector{Chromosome},Vector{LinCircuit}}=Chromosome[] )
  if length(pop) == 0
    pop = random_population( p, funcs, popsize, use_cgp=use_cgp )
  end
  fitness_vector = Float64[]  # Establish scope
  first_opt_gen = 0
  first_majority_gen = 0
  for gen = 1:ngens
    mutate_population!( pop, mutrate, funcs )
    fitness_vector = map( ch->fitness_funct( ch, target ), pop )
    #println("A gen: ",gen,"  fitness_vector: ",fitness_vector)
    for i = 1:popsize
      fit = fitness_funct( pop[i], target )
      pop[i].fitness = fitness_vector[i]
    end
    maxfit = findmax(fitness_vector)[1]
    first_opt_gen = (first_opt_gen==0 && maxfit == 1.0) ? gen : first_opt_gen
    first_majority_gen = (first_majority_gen==0) && (length( findall( x->x==1.0, fitness_vector )) >= ceil(popsize/2)) ? gen : first_majority_gen
    propsel!( pop, fitness_vector, maxfit=maxfit )
    fitness_vector = map( ch->fitness_funct( ch, target ), pop )
    #println("B gen: ",gen,"  fitness_vector: ",fitness_vector)
  end
  meanfit = mean( fitness_vector )
  maxfit = findmax(fitness_vector)[1]
  #pop
  ( meanfit, maxfit, first_opt_gen, first_majority_gen )
end

function fitness_funct( c::Chromosome, target::Goal )
  p = c.params
  1.0-hamming_distance( target, output_values(c), p.numinputs ) 
end

function fitness_funct( c::Chromosome, target_list::GoalList )
  p = c.params
  maximum( 1.0-hamming_distance( g, output_values(c), p.numinputs ) for g in target_list )
end

# mutrate is the probability that a population member will be point-mutated in one generation.
# This does allow for the possibility of more than one mutation to a population member in a generation
function mutate_population!( pop::Union{Vector{Chromosome},Vector{LinCircuit}}, mutrate::Float64, funcs::Vector{Func} )
  use_cgp = (typeof(pop) == Vector{Chromosome})
  p = pop[1].params
  Pois = Poisson( length(pop)*mutrate )
  nmuts = rand(Pois)
  #println("nmuts: ",nmuts)
  for i = 1:nmuts
    j = rand(1:length(pop))
    old_ov = output_values(pop[j])
    if use_cgp
      mutate_chromosome!( pop[j], funcs )
    else
      mutate_circuit!( pop[j], funcs )
    end
    new_ov = output_values(pop[j])
    #println("mutated pop[",j,"]  old_ov: ",old_ov,"  new_ov: ",new_ov)
  end
  #println(map(output_values,pop))
end

function random_population( p::Parameters, funcs::Vector{Func}, popsize::Int64; use_cgp::Bool=true )
  if use_cgp
    return [ random_chromosome( p, funcs) for _=1:popsize ]
  else
    return [ rand_lcircuit( p, funcs) for _=1:popsize ]
  end
end
