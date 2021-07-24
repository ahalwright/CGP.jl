function run_run_pop_evolvability( nreps::Int64, p::Parameters, popsize_rng::CGP.IntRange, gl::GoalList, ngens_popsize::Int64, mutrate_rng::CGP.FloatRange;
    csvfile::String="", prdebug::Bool=false )
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
  results=pmap(x->run_pop_evolvability(nreps,p,x[1],gl,x[2],x[3]),run_tuples)
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
      println(f,"# nreps: ",nreps)
      CSV.write( f, sumdf, append=true, writeheader=true )
    end
  end
  sumdf
end
      
function run_pop_evolvability( nreps::Int64, p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, mutrate::Float64=1.0;
    csvfile::String="", prdebug::Bool=false )
  println("mutrate: ",mutrate)
  avgdf = DataFrame()
  avgdf.gen = collect(1:ngens)
  avgdf.popsize = fill(popsize,ngens)
  avgdf.mutrate = fill(mutrate,ngens)
  avgdf.max_fitness = zeros(Float64,ngens)
  avgdf.evolvability = zeros(Float64,ngens)
  avgdf.robustness = zeros(Float64,ngens)
  #dflist = pmap(x->pop_evolvability( p, popsize, gl, ngens, mutrate, prdebug=prdebug ), collect(1:nreps))
  dflist = map(x->pop_evolvability( p, popsize, gl, ngens, mutrate, prdebug=prdebug ), collect(1:nreps))
  for df in dflist
    # df = pop_evolvability( p, popsize, gl, ngens, mutrate, prdebug=prdebug )
    avgdf.max_fitness += df.max_fitness
    avgdf.evolvability += df.evolvability
    avgdf.robustness += df.robustness
  end
  avgdf.max_fitness /= nreps
  avgdf.evolvability /= nreps
  avgdf.robustness /= nreps
  avgdf[end,:]
end

function pop_evolvability( p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, mutrate::Float64=1.0; 
    csvfile::String="", prdebug::Bool=false )
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
    c.fitness = fitness_funct( p, c, gl )
    prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
    fitness_vector[i] = pop[i].fitness 
    maxfit = fitness_vector[i] > maxfit ? fitness_vector[i] : maxfit
  end   
  ( evolvability, robustness ) = pop_evolvability_robustness( pop )
  push!( df.gen, 1 )
  push!( df.max_fitness, maxfit )
  push!( df.evolvability, evolvability )
  push!( df.robustness, robustness )
  for gen = 2:ngens
    for i in 1:popsize
      c = pop[i]
      if rand() <= mutrate
        pop[i] = c = mutate_chromosome!(deepcopy(c),funcs)[1]
        c.fitness = fitness_funct( p, c, gl )
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
  println("gl: ",gl,"  max_fit: ",maxfit)
  df
end

function fitness_funct( p::Parameters, c::Chromosome, gl::GoalList )
  maximum( 1.0-hamming_distance( g, output_values(c), p.numinputs ) for g in gl )
end
