# Test whether degeneracy increases with numgates (numinteriors).
# Test whether degeneracy increases with sporadic mutation of the goal.
export recover_phenotype

function test_degen( numinputs::Int64, numoutputs::Int64, numints_rng::IntRange,
    numtrials::Int64, maxsteps::Int64 )
  df = DataFrame()
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numints = Int64[]
  df.levsback = Int64[]
  df.degen = Float64[]
  for numints in numints_rng
    levsback = Int(ceil(numints/2))
    sum_degen = 0.0
    for t = 1:numtrials
      p = Parameters( numinputs, numoutputs, numints, levsback )
      c = random_chromosome(p)
      g = randgoal(numinputs,numoutputs)
      new_c = neutral_evolution( c, g, maxsteps )
      sum_degen += degeneracy( c )
    end
    row = (numinputs,numoutputs,numints,levsback,sum_degen/numtrials)
    push!(df,row)
  end
  df
end

# Helper function for test_degen()
function run_test_degen(numinputs::Int64, numoutputs::Int64, numints::Int64, maxsteps::Int64 )
  p = Parameters( numinputs, numoutputs, numints, levsback )
  c = random_chromosome(p)
  g = randgoal(numinputs,numoutputs)
  new_c = neutral_evolution( c, g, maxsteps )
  sum_degen += degeneracy( c )
end

# Starts with a chromsome c.
# For each trial, mutate the chromosome, and use neutral evolution to try evolve the mutated chromosome
#   to output the phenotype of c.  Up to maxtries attempts at neutral evolution are made
# Results the are the sum of the evolutionary steps made by neutral evolution
#   and the sum of the number of tries at neutral evolution. 
function recover_phenotype( c::Chromosome, maxsteps::Int64, maxtrials::Int64, maxtries::Int64 )
  p = c.params
  funcs = default_funcs(p.numinputs)
  g = output_values(c)
  sumsteps = 0
  sumtries = 0
  for t = 1:maxtrials
    (mc,active) = mutate_chromosome!(deepcopy(c),funcs)
    (nc,steps) = neutral_evolution( mc, g, maxsteps )
    sumsteps += steps
    tries = 1
    while tries < maxtries && steps == maxsteps
      (nc,steps) = neutral_evolution( mc, g, maxsteps )
      sumsteps += steps
      tries += 1
    end 
    if tries == maxtries
      println("tries: ",tries,"  maxtries: ",maxtries)
    end
    sumtries += tries
  end
  return (sumsteps,sumtries)
end

#=
# Moved to Evolve.jl.
function lambda_evolution( c::Chromosome, g::Goal, maxsteps::Integer, mutrate::Float64 )
  p = c.params
  p.mutrate = mutrate
  funcs = default_funcs(p.numinputs)
  poisson_lambda = p.mutrate*p.numinteriors
  X = Poisson( poisson_lambda )
  step = 0
  ov = output_values( c)
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  #print("ov: ",ov,"  cur_dist: ",current_distance,"  "); print_circuit(c)
  while step < maxsteps && ov != g
    chrome_list = Chromosome[]
    dist_list = Float64[]
    for i = 1:p.lambda
      step += 1
      num_mutations =  rand(X)
      new_c = deepcopy(c) 
      for m = 1:num_mutations
        mutate_chromosome!( new_c, funcs )
      end
      new_ov = output_values( new_c)
      new_dist = hamming_distance( new_ov, g, c.params.numinputs )
      #print("new_ov: ",new_ov,"  new_dist: ",new_dist,"  "); print_circuit(new_c)
      push!(chrome_list,new_c)
      push!(dist_list,new_dist)
    end  
    (best_dist,ind) = findmin( dist_list )
    if best_dist <= current_distance
      c = chrome_list[ind]
      ov = output_values( c )
      current_distance = hamming_distance( ov, g, c.params.numinputs )
      #print("ov: ",ov,"  cur_dist: ",current_distance,"  "); print_circuit(c)
    end
  end # while
  if step == maxsteps
    println("lambda evolution failed with ",step," steps for goal: ",g)
    return (c, step)
  else
    println("lambda evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end        
=#

function compare_lambda_neutral( p::Parameters, g::Goal, trials::Int64, maxsteps::Int64, mutrate_range::CGP.FloatRange )
  @assert length(collect(mutrate_range)) == 3
  df = DataFrame()
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numinteriorss = Int64[]
  df.numlevelsback = Int64[]
  df.mutrate = Float64[]
  df.neutral_fails = Int64[]
  df.neutral_steps= Int64[]
  df.lamda1_fails = Int64[]
  df.lamda1_steps = Int64[]
  df.lamda2_fails = Int64[]
  df.lamda2_steps = Int64[]
  df.lamda3_fails = Int64[]
  df.lamda3_steps = Int64[]
  for t = 1:numtrials
    c = random_chromosome(p)
    (neut_c,neut_steps) = neutral_evolution( deepcopy(c), g, maxsteps )
    mr = collect(mutrate_range)
    (lambda_c1,lambda_steps1) = lambda_evolution( deepcopy(c), g, maxsteps, mr[1] )
    (lambda_c2,lambda_steps2) = lambda_evolution( deepcopy(c), g, maxsteps, mr[2] )
    (lambda_c3,lambda_steps3) = lambda_evolution( deepcopy(c), g, maxsteps, mr[3] )
    df_row = (
      p.numinputs,
      p.numoutputs
    )
  end
end
      

