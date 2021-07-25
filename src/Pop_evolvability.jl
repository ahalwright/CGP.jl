function run_run_pop_evolvability( nreps::Int64, p::Parameters, popsize_rng::CGP.IntRange, gl::GoalList, ngens_popsize::Int64, mutrate_rng::CGP.FloatRange;
    csvfile::String="" )
  sumdf = DataFrame()
  sumdf.gen = Int64[]
  sumdf.popsize = Int64[]
  sumdf.mutrate = Float64[]
  sumdf.max_fitness = Float64[]
  sumdf.evolvability = Float64[]
  sumdf.robustness= Float64[]
  run_tuples = Tuple[]
  for popsize in popsize_rng
    ngens = Int(round(ngens_popsize/popsize))
    for mutrate in mutrate_rng
      push!(run_tuples, (popsize,ngens,mutrate))
    end
  end
  results=pmap(x->run_pop_evolvability(nreps,p,x[1],gl,x[2],x[3],use_pmap=false),run_tuples) 
  for res in results
    push!(sumdf,res)
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# goallist: ",gl)
      println(f,"# nreps: ",nreps)
      CSV.write( f, sumdf, append=true, writeheader=true )
    end
  end
  sumdf
end
      
# Averages results over reps
# use_pmap must be false when called by run_run_pop_evolvability
# If length(gl) == 0, generate goallist within each call to pop_evolvability
function run_pop_evolvability( nreps::Int64, p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, 
    start_phmut_gen::Int64, phmut_interval::Int64,  mutrate::Float64=1.0; numgoals::Int64=1,
    csvfile::String="", use_pmap::Bool=false, all_gens::Bool=false )
  #@assert p.numoutputs == length(gl[1])
  println("mutrate: ",mutrate,"  gl: ",gl,"  length(gl): ",length(gl))
  avgdf = DataFrame()
  avgdf.gen = collect(1:ngens)
  avgdf.popsize = fill(popsize,ngens)
  avgdf.mutrate = fill(mutrate,ngens)
  avgdf.max_fitness = zeros(Float64,ngens)
  avgdf.evolvability = zeros(Float64,ngens)
  avgdf.robustness = zeros(Float64,ngens)
  if use_pmap
    dflist = pmap(x->pop_evolvability( p, popsize, gl, ngens, start_phmut_gen, phmut_interval, mutrate,numgoals=numgoals), collect(1:nreps))
  else
    dflist = map(x->pop_evolvability( p, popsize, gl, ngens, start_phmut_gen, phmut_interval, mutrate,numgoals=numgoals), collect(1:nreps))
  end
  for df in dflist
    avgdf.max_fitness += df.max_fitness
    avgdf.evolvability += df.evolvability
    avgdf.robustness += df.robustness
  end
  avgdf.max_fitness /= nreps
  avgdf.evolvability /= nreps
  avgdf.robustness /= nreps
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# goallist: ",gl)
      println(f,"# nreps: ",nreps)
      CSV.write( f, avgdf, append=true, writeheader=true )
    end
  end
  if all_gens
    avgdf
  else
    avgdf[end,:]
  end
end

# If length(gl)==0, generate goallist within each call to pop_evolvability
function pop_evolvability( p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, 
    start_phmut_gen::Int64, phmut_interval::Int64,  mutrate::Float64=1.0; 
    csvfile::String="", numgoals::Int64=1 )
  if length(gl) > 0
    glmut = deepcopy(gl)
  else
    glmut = randgoallist( numgoals, p )
  end
  println("glmut: ",glmut)
  df = DataFrame()
  df.gen = Int64[]
  df.max_fitness = Float64[]
  df.evolvability = Float64[]
  df.robustness = Float64[]
  funcs = default_funcs( p.numinputs )
  pop = [ random_chromosome(p,funcs) for i = 1: popsize ]
  count_mutations = 0
  fitness_vector = zeros(Float64,popsize)
  maxfit = 0.0
  for i in 1:popsize
    c = pop[i]
    c.fitness = fitness_funct( p, c, glmut )
    fitness_vector[i] = pop[i].fitness 
    maxfit = fitness_vector[i] > maxfit ? fitness_vector[i] : maxfit
  end   
  ( evolvability, robustness ) = pop_evolvability_robustness( pop )
  push!( df.gen, 1 )
  push!( df.max_fitness, maxfit )
  push!( df.evolvability, evolvability )
  push!( df.robustness, robustness )
  for gen = 2:ngens
    if gen >= start_phmut_gen && gen % phmut_interval == 0
      phmut!( glmut, p )
      println("phmut gen: ",gen,"  glmut: ",glmut)
    end
    for i in 1:popsize
      c = pop[i]
      if rand() <= mutrate
        pop[i] = c = mutate_chromosome!(deepcopy(c),funcs)[1]
        c.fitness = fitness_funct( p, c, glmut )
        #println("mutation count: ",count_mutations,"  phenotype: ",@sprintf("0x%x",(output_values(c)[1])),"  fitness: ",c.fitness)
        count_mutations += 1
      end
    end   
    fitness_vector = ([ pop[i].fitness for i = 1:popsize ])
    maxfit = findmax(fitness_vector)[1]
    propsel!( pop, fitness_vector, maxfit=maxfit )
    ( evolvability, robustness ) = pop_evolvability_robustness( pop )
    push!( df.gen, gen )
    push!( df.max_fitness, maxfit )
    push!( df.evolvability, evolvability )
    push!( df.robustness, robustness )
  end
  println("glmut: ",glmut,"  max_fit: ",maxfit)
  df
end

function phmut!( glmut::GoalList, p::Parameters ) 
  gl_component = length(glmut) > 1 ? rand(1:length(glmut)) : 1
  goal_index = length(glmut[1]) > 1 ? rand(1:length(glmut[1])) : 1
  max_shift = 2^p.numinputs - 1
  shift = rand(0:max_shift)
  println("gl_component: ",gl_component,"  goal_index: ",goal_index,"  shift: ",shift)
  mut_mask = MyInt(1) << shift
  #@printf("mut_mask: 0x%x\n",mut_mask)
  # A more elegant version of the next line is:  
  # glmut[gl_component][goal_index] ⊻= mut_mask
  glmut[gl_component][goal_index] = xor(glmut[gl_component][goal_index], mut_mask)
  return glmut
end


function fitness_funct( p::Parameters, c::Chromosome, gl::GoalList )
  maximum( 1.0-hamming_distance( g, output_values(c), p.numinputs ) for g in gl )
end
