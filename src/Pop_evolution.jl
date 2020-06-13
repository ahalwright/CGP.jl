# Implements population based evoluiont with runs of individual evolution in each population generation
# Populations are lists of chromosomes
export pop_evolve

function pop_evolve( params::Parameters, popsize::Int64, max_pop_steps::Int64, max_indiv_steps::Int64,
    goallist::GoalList, funcs::Vector{Func}; uniform_start::Bool=false, hamming_sel::Bool=true )
  if uniform_start
    c = random_chromosome(params,funcs)
    pop = fill( c, popsize )
  else
    pop = fill( random_chromosome(params,funcs), popsize )
  end
  max_gen_fitness = zeros(Float64,max_pop_steps)
  for g = 1:max_pop_steps
    println("starting generation ",g,"  fitnesses: ",[pop[i].fitness for i =1:popsize])
    for i = 1:popsize
      result = mut_evolve( pop[i], goallist, funcs, max_indiv_steps, hamming_sel=hamming_sel ) 
      if length(result) == 2
        println("returning (orig_c,c)")
        return result
      else
        (pop[i],step,worse,same,better,output,matched_goals_list) = result
      end
      #(pop[i],step,worse,same,better,output,matched_goals_list) = mut_evolve( pop[i], goallist, funcs, max_indiv_steps, hamming_sel=hamming_sel )
      println("(g,i): ",(g,i),"  pop[i].fitness: ",pop[i].fitness)
    end
    # Do proportional selection
    fitness_vector = [ pop[i].fitness for i = 1:popsize ]
    println("gen: ",g,"  max fit before selection: ",findmax(fitness_vector)[1])
    (max_gen_fitness[g],maxfit_index) = findmax(fitness_vector)
    if max_gen_fitness[g] == params.numoutputs
      break  # found optimum so break out of loop
    end
    propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] ) 
    fitness_vector = [ pop[i].fitness for i = 1:popsize ]
    (max_gen_fitness[g],maxfit_index) = findmax(fitness_vector)
    println("gen: ",g,"  max fit after selection: ",max_gen_fitness[g],"  at index: ",maxfit_index)
  end
  max_fit_all_generations = findmax(max_gen_fitness)[1]
  println("max fit over generations: ",max_fit_all_generations)
  return (max_fit_all_generations,pop)
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

