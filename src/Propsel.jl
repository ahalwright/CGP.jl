#=
Functions to implement proportional selection.
Returns a new population where the individuals are selected in proportion to their fitness.
=#

export  propsel, propsel!
@doc """ function propsel( pop::Population, fitness::Function )
Apply proportional selection to Population pop using fitness, 
and return the result.  
"""
function propsel( pop::Population, dfe::Function, fitness_table::Dict{Int64,Float64} )
  new_pop = deepcopy(pop)
  propsel!( new_pop, dfe, fitness_table )
  new_pop
end

@doc """ function propsel( pop::Population, fitness::Function )
Apply proportional selection to Population pop using fitness to generate
and return a new selected population of size n.
"""
function propsel( pop::Population, n::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
  new_pop = zeros(Int64,n)
  fmax = 0.0
  for p in pop
    fitp = dfe_fitness( p, dfe, fitness_table )
    if fitp > fmax
      fmax = fitp
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    return
  end

  N = length(pop)
  k = 0
  while k < n
    i = rand(1:N)
    w = dfe_fitness(pop[i],dfe, fitness_table ) / fmax
    if rand() < w
      new_pop[k + 1] = pop[i]
      k += 1
    end
  end
  new_pop
end

#=
@doc """function propsel!(p::Population, fitness::Vector{Float64} )
Conduct proportional selection in-place.
"""
function propsel!( pop::Population, fitness::Vector{Float64}  )
  fmax = maximum(fitness)
  if fmax == 0
    # all elements have fitness zero
    return
  end

  n = length(pop)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = fitness[pop[i]] / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end
=#

@doc """function propsel!(p::Population, fitness::Function, dfe::Function )
Conduct proportional selection in-place.
"""
function propsel!( pop::Population, dfe::Function, fitness_table::Dict{Int64,Float64} )
  fmax = 0.0
  for p in pop
    fitp = dfe_fitness( p, dfe, fitness_table )
    if fitp > fmax
      fmax = fitp
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    return
  end

  n = length(pop)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = dfe_fitness(pop[i], dfe, fitness_table ) / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end

@doc """ function propsel_alt()
  Roulette wheel version of proportial selection
  Less efficient in almost all circumstances than the above version.
  Written as check of correctness for the above version
"""
function propsel_alt( pop::Population, dfe::Function)
  fitness_table = Dict{Int64,Float64}()
  N = length(pop)
  new_pop = zeros(Int64,N)
  fitdict = Dict{Int64,Float64}()
  popctr = pop_counter( pop )
  for p in pop
    fit = dfe_fitness(p, dfe, fitness_table )  # set fitness of p to be dfe(p).
  end
  sum_fitness = 0.0
  for p in keys(popctr)
    dfit = dfe_fitness( p, dfe, fitness_table )*popctr[p]
    sum_fitness += dfit
  end
  for p in keys(popctr)
    fitdict[p] = dfe_fitness( p, dfe, fitness_table )*popctr[p]/sum_fitness
  end
  for i = 1:N
    r = rand()
    sumfit = 0.0
    for p in keys(fitdict)
      sumfit += fitdict[p]
      if r < sumfit
        new_pop[i] = p
        break
      end
    end
  end
  new_pop
end

