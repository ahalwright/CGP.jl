# Population-based evolution to try to approximately replicate Labar & Adami 2016.
# The Julia type PredType is defined in Chromosome.jl
using Main.CGP
using Statistics

# mutrate is the probability that a chromosome (circuit) will be mutation in a generation
function run_both_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, goallist::GoalList, 
    ngens::Int64, numpops::Int64, mutrate::Float64=1.0;
    csvfile::String="", uniform_start::Bool=true, prdebug::Bool=false )
  df=DataFrame()
  df.goal = Goal[]
  df.ngens = Int64[]
  df.numpops = Int64[]
  df.popsize = Int64[]
  df.mutrate = Float64[]
  df.evolvability = Int64[]
  df.maxfit = Float64[]
  df.maxfit_gen = Int64[]
  df.fit_decreases = Int64[]
  df_list = pmap(goal->run_both_one_goal(nreps,p,popsize,goal,ngens,numpops,mutrate,uniform_start=uniform_start),goallist)
  #df_list = map(goal->run_both_one_goal(nreps,p,popsize,goal,ngens,numpops,mutrate,uniform_start=uniform_start),goallist)
  df = vcat( df_list... )
  #=
  for goal in goallist
    println("goal: ",goal)
    df = run_pop_evolve(nreps,p,popsize,goal,ngens,df=df,uniform_start=uniform_start)
    mpopsize = Int64(round(popsize/numpops))   # popsize should be divisible by numpops
    df = run_multiple_pop_evolve(nreps,p,mpopsize,goal,ngens,numpops,mutrate,df=df,uniform_start=uniform_start)
  end
  =#
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# ngens: ",ngens)
      println(f,"# nreps: ",nreps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function run_both_one_goal( nreps::Int64, p::Parameters, popsize::Int64, goal::Goal, ngens::Int64, numpops::Int64,
    mutrate::Float64=1.0;
    uniform_start::Bool=true, prdebug::Bool=false )
  println("run both one goal: ",goal)
  df = DataFrame()
  df.goal = Goal[]
  df.ngens = Int64[]
  df.numpops = Int64[]
  df.popsize = Int64[]
  df.mutrate = Float64[]
  df.evolvability = Int64[]
  df.maxfit = Float64[]
  df.maxfit_gen = Int64[]
  df.fit_decreases = Int64[]
  println("goal: ",goal)
  df = run_pop_evolve(nreps,p,popsize,goal,ngens,mutrate,df=df,uniform_start=uniform_start)
  #println("run_both_one_goal A size(df): ",size(df))
  mpopsize = Int64(round(popsize/numpops))   # popsize should be divisible by numpops
  df = run_multiple_pop_evolve(nreps,p,mpopsize,goal,ngens,numpops,mutrate,df=df,uniform_start=uniform_start)  
  #println("run_both_one_goal B size(df): ",size(df))
  df
end

function run_multiple_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, g::Goal, ngens::Int64, numpops::Int64,
    mutrate::Float64=1.0;
    df::DataFrame=DataFrame(), csvfile::String="", uniform_start::Bool=true, prdebug::Bool=false )
  println("run_multiple_pop_evolve goal: ",g)
  if size(df)[1] == 0
    df.goal = Goal[]
    df.ngens = Int64[]
    df.numpops = Int64[]
    df.popsize = Int64[]
    df.mutrate = Float64[]
    df.evolvability = Int64[]
    df.maxfit = Float64[]
    df.maxfit_gen = Int64[]
    df.fit_decreases = Int64[]
  end
  for rep = 1:nreps
    (mutrate,mphenotypes,mpredecessors,mfitnesses,mcomplexities,mnumgates)=multiple_pop_evolve(p,popsize,g,ngens,numpops,mutrate,uniform_start=uniform_start)
    #println("run_multiple_pop_evolve size(mphenotypes): ",size(mphenotypes),"  size(mphenotypes[1]): ",size(mphenotypes[1]))
    push!(df.goal,g)
    push!(df.ngens,ngens)
    push!(df.numpops,numpops)
    push!(df.popsize,popsize)
    push!(df.mutrate,mutrate)
    push!(df.evolvability,final_evolvability(mphenotypes))
    findmax_mfh=findmax(max_fitness_history(mfitnesses))
    push!(df.maxfit,findmax_mfh[1][1])
    push!(df.maxfit_gen,findmax_mfh[2])
    push!(df.fit_decreases,preds_fit_decreases( mpredecessors, mphenotypes, mfitnesses )[1] )
  end
  df
end
 
function run_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, g::Goal, ngens::Int64, mutrate::Float64=1.0;
    df::DataFrame=DataFrame(), csvfile::String="", uniform_start::Bool=true, prdebug::Bool=false )
  println("run_pop_evolve goal: ",g)
  if size(df)[1] == 0
    df.goal = Goal[]
    df.ngens = Int64[]
    df.numpops = Int64[]
    df.popsize = Int64[]
    df.mutrate = Float64[]
    df.evolvability = Int64[]
    df.maxfit = Float64[]
    df.maxfit_gen = Int64[]
    df.fit_decreases = Int64[]
  end
  for rep = 1:nreps
    (pop,phenotypes,predecessors,fitnesses,complexities,numgates)=
        pop_evolve(p,popsize,g,ngens,mutrate,uniform_start=uniform_start)
    #println("run_pop_evolve size(phenotypes): ",size(phenotypes))
    push!(df.goal,g)
    push!(df.ngens,ngens)
    push!(df.numpops,1)
    push!(df.popsize,popsize)
    push!(df.mutrate,mutrate)
    push!(df.evolvability,final_evolvability(phenotypes))
    findmax_mfh=findmax(max_fitness_history(fitnesses))
    push!(df.maxfit,findmax_mfh[1])
    push!(df.maxfit_gen,findmax_mfh[2])
    push!(df.fit_decreases,preds_fit_decreases( predecessors, phenotypes, fitnesses )[1] )
    if nreps==1
      return (df,pop,phenotypes,predecessors,fitnesses,complexities,numgates)
    end
  end
  df
end
    
function multiple_pop_evolve( p::Parameters, popsize::Int64, g::Goal, ngens::Int64, numpops::Int64, mutrate::Float64=1.0;
    uniform_start::Bool=true, prdebug::Bool=false )
  mphenotypes = fill(zeros(MyInt,0,0),numpops)
  mpredecessors = fill(zeros(PredType,0,0),numpops)
  mfitnesses = fill(zeros(Float64,0,0),numpops)  
  mcomplexities = fill(zeros(Float64,0,0),numpops)
  mnumgates = fill(zeros(Int64,0,0),numpops)
  for pind = 1:numpops
    #println("pind: ",pind)
    (pop,mphenotypes[pind],mpredecessors[pind],mfitnesses[pind],mcomplexities[pind],mnumgates[pind]) = 
        pop_evolve(p,popsize,g,ngens,mutrate,uniform_start=false)
  end
  (mutrate,mphenotypes,mpredecessors,mfitnesses,mcomplexities,mnumgates)
end

function pop_evolve( p::Parameters, popsize::Int64, g::Goal, ngens::Int64, mutrate::Float64=1.0; 
    uniform_start::Bool=true, prdebug::Bool=false )
  #println("pop_evolve: goal: ",g,"  psize: ",popsize,"  ngens: ",ngens,"  mutrate: ",mutrate)
  count_mutations = 0
  funcs = default_funcs( p.numinputs )
  phenotypes = zeros(MyInt,ngens,popsize)
  predecessors = zeros(PredType,ngens,popsize)
  fitnesses = zeros(Float64,ngens,popsize)
  complexities = zeros(Float64,ngens,popsize)
  numgates = zeros(Int64,ngens,popsize)
  if uniform_start
    c = random_chromosome(p,funcs,ident=PredType(1))
    pop = fill( c, popsize )
  else
    pop = [ random_chromosome(p,funcs,ident=PredType(i)) for i = 1: popsize ]
  end
  for i in 1:popsize
    c = pop[i]
    c.robustness = i
    c.fitness = 1.0-hamming_distance( g, output_values(c), p.numinputs )
    prdebug ? print(c.robustness,"  ") : nothing
    prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
  end
  for gen = 1:ngens
    #fitness_vector = rescale_fitnesses([ pop[i].fitness for i = 1:popsize ])
    fitness_vector = ([ pop[i].fitness for i = 1:popsize ])
    prdebug ? println("fit_vect: ",fitness_vector) : nothing
    prdebug ? println("gen: ",gen,"  max fit: ",maximum(fitness_vector)) : nothing
    propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] )
    prdebug ? println("after propsel gen: ",gen) : nothing
    for c in pop
      c.fitness = hamming_distance( g, output_values(c), p.numinputs )
      prdebug ? print(c.robustness,"  ") : nothing
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
    end
    #=
    for i in 1:length(pop)
      phenotypes[gen,i] = output_values(pop[i])[1]
      predecessors[gen,i] = pop[i].robustness
      fitnesses[gen,i] = pop[i].fitness
      complexities[gen,i] = complexity5(pop[i])
      numgates[gen,i] = pop[i].params.numinteriors
    end
    if prdebug
      println("phenos:  ",phenotypes[gen,:])
      println("preds:   ",predecessors[gen,:])
      println("fitness: ",fitnesses[gen,:])
      println("complex: ",complexities[gen,:])
      #println("numgates:",numgates[gen,:])
    end
    =#
    sav_pop = deepcopy( pop )
    for i in 1:popsize
      c = pop[i]
      if rand() <= mutrate
        #print("befor mut: ")
        #print_circuit(c)
        pop[i] = c = mutate_chromosome!(deepcopy(c),funcs)[1] 
        #print("after mut: ")
        #print_circuit(c)
        #print("newc:      ")
        #print_circuit(new_c)
        c.fitness = 1.0-hamming_distance( g, output_values(c), p.numinputs )
        #println("mutation count: ",count_mutations,"  phenotype: ",@sprintf("0x%x",(output_values(c)[1])),"  fitness: ",c.fitness)
        count_mutations += 1
      end
      c.robustness = i  # Robustness is not robustness; instead it is an identifier of the index of circut in the population
    end
    
    for i in 1:length(pop)
      phenotypes[gen,i] = output_values(pop[i])[1]
      predecessors[gen,i] = pop[i].robustness
      fitnesses[gen,i] = pop[i].fitness = 1.0-hamming_distance( g, output_values(pop[i]), p.numinputs )
      @assert pop[i].fitness == 1.0-hamming_distance( g, output_values(pop[i]), p.numinputs )  
      complexities[gen,i] = complexity5(pop[i])
      numgates[gen,i] = pop[i].params.numinteriors
    end
    if prdebug
      println("phenos:  ",phenotypes[gen,:])
      println("preds:   ",predecessors[gen,:])
      println("fitness: ",fitnesses[gen,:])
      println("complex: ",complexities[gen,:])
      #println("numgates:",numgates[gen,:])
    end
    
    prdebug ? println("after mutate gen: ",gen) : nothing
    #=
    for i in 1:popsize
      c = pop[i]
      c.robustness = i  # Robustness is not robustness; instead it is an identifier of the index of circut in the population
      c.fitness = 1.0-hamming_distance( g, output_values(c), p.numinputs )
      prdebug ? print(c.robustness,"  ") : nothing
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
    end
    =#
  end
  println("end pop evolve: len(pop): ",length(pop),"  size: ",size(phenotypes))
  (pop,phenotypes,predecessors,fitnesses,complexities,numgates)
end

function rescale_fitnesses( fit_vect::Vector{Float64} )
  fit_min = minimum( fit_vect )
  #fit_min = quantile( fit_vect, 0.20 )
  fit_max = maximum( fit_vect )
  frange = fit_max - fit_min
  if frange > 0.0
    return [ (f-fit_min)/frange for f in fit_vect ]
  else
    return fit_vect
  end
end

function pheno_history( phenotypes::Array{MyInt,2} )
  (num_gens, popsize) = size(phenotypes)
  pheno_history = fill(Set(MyInt[]),num_gens)
  pheno_set = Set(MyInt[])
  for g = 1:num_gens
    union!(pheno_set,Set(phenotypes[g,:]))
    pheno_history[g] = deepcopy(pheno_set)
  end
  pheno_history
end

# Result should agree with evolvability()
function pheno_history( mphenotypes::Array{Array{MyInt,2},1} )
  npops = length(mphenotypes)
  (num_gens, popsize) = size(mphenotypes[1])
  #phenohistories = [pheno_history( mphenotypes[i] ) for i = npops ]  # fails due to a bug in Julia
  phenohistories =  Array{Set{UInt16},1}[]
  for i = 1:npops
    push!(phenohistories,pheno_history( mphenotypes[i]))
  end
  phenoHistory = fill(Set(MyInt[]),num_gens)
  for g = 1:num_gens
    phenoHistory[g] = phenohistories[1][g]
    for i = 1:npops
      phenoHistory[g] = union!(phenoHistory[g],phenohistories[i][g])
    end
  end
  pheno_set = Set(MyInt[])
  for g = 1:num_gens
    union!(pheno_set,Set(phenoHistory[g]))
    phenoHistory[g] = deepcopy(pheno_set)
  end
  phenoHistory
end

function final_pheno_history( mphenotypes::Array{Array{MyInt,2},1} )
  length(pheno_history(mphenotypes)[end])
end

function evolvability( phenotypes::Array{MyInt,2} )
  (num_gens, popsize) = size(phenotypes)
  evo_history = zeros(Int64,num_gens)
  pheno_set = Set(MyInt[])
  for i = 1:num_gens
    union!(pheno_set,Set(phenotypes[i,:]))
    evo_history[i] = length(pheno_set)
  end
  pheno_set
  #evo_history
end

function final_evolvability( phenotypes::Array{MyInt,2} )
  length(evolvability(phenotypes))
end

# Result should agree with pheno_history()
function evolvability( mphenotypes::Array{Array{MyInt,2},1} )
  npops = length(mphenotypes)
  (num_gens, popsize) = size(mphenotypes[1])
  # Take the union of the phenotypes of each population for each generation g
  ph_gen = [union( map( union, [ mphenotypes[j][g,:] for j = 1:npops] )... ) for g = 1:num_gens]
  pheno_history = fill(Set(MyInt[]),num_gens)
  pheno_set = Set(MyInt[])
  for g = 1:num_gens
    union!(pheno_set,Set(ph_gen[g]))
    pheno_history[g] = deepcopy(pheno_set)
  end
  pheno_history
end

function final_evolvability( mphenotypes::Array{Array{MyInt,2},1} )
  length(evolvability(mphenotypes)[end])
end

function max_fitness_history( fitnesses::Array{Float64,2} )
  (num_gens, popsize) = size(fitnesses)
  history = [maximum( fitnesses[g,:]) for g = 1:num_gens ]
end

function max_fitness_history( mfitnesses::Array{Array{Float64,2},1} )
  numpops = length(mfitnesses)
  (num_gens, popsize) = size(mfitnesses[1])
  histories = [[ maximum(mfitnesses[i][g,:]) for g = 1:num_gens ] for i = 1:numpops ]
  history = [findmax( [histories[i][g] for i = 1:numpops] ) for g = 1:num_gens ]
end

function mean_fitness_history( fitnesses::Array{Float64,2} )
  (num_gens, popsize) = size(fitnesses)
  history = [mean( fitnesses[g,:]) for g = 1:num_gens ]
end
  
function max_complexity_history( complexities::Array{Float64,2} )
  (num_gens, popsize) = size(complexities)
  history = [maximum( complexities[g,:]) for g = 1:num_gens ]
end

function mean_complexity_history( complexities::Array{Float64,2} )
  (num_gens, popsize) = size(complexities)
  history = [mean( complexities[g,:]) for g = 1:num_gens ]
end
  
# Return the list of unique predecessors of element index of populations 1 to gen.
# Usually one will call with gen = maxgenerations
function preds( predecessors::Array{Int64,2}, gen::Int64, index::Int64 )
  indices = zeros(Int64,gen)
  while gen > 0
    indices[gen] = index = predecessors[gen,index]
    #println("gen: ",gen,"  index: ",index)
    gen -= 1
  end
  indices
end

# Returns a list of phenotype-fitness pairs of the predecessors of a random maximum fitness element
#   of the first population to reach the maximum fitness.
function preds_phenos_fits( predecessors::Array{Int64,2}, phenotypes::Array{MyInt,2}, fitnesses::Array{Float64,2} ) 
  #println("preds_phenos_fits: size(phenotypes): ",size(phenotypes))
  (maxfit,maxfit_gen) = findmax_mfh=findmax(max_fitness_history(fitnesses))
  column = rand(findmaxall( fitnesses[maxfit_gen,:])[2])  # A random columm with maximum fitness
  prds = preds( predecessors, maxfit_gen, column )                         
  [(phenotypes[i,prds[i]],fitnesses[i,prds[i]]) for i = 1:maxfit_gen ]
end

function preds_phenos_fits( mpredecessors::Array{Array{Int64,2},1}, mphenotypes::Array{Array{MyInt,2},1}, 
    mfitnesses::Array{Array{Float64,2},1} ) 
  #println("preds_phenos_fits  size(mphenotypes[1]): ",size(mphenotypes[1]))
  ((maxfit,maxfit_pop),maxfit_gen) = findmax_mfh=findmax(max_fitness_history(mfitnesses))
  #println("preds_phenos_fits  size(mfitnesses[maxfit_pop]): ",size(mfitnesses[maxfit_pop]))
  column = rand(findmaxall( mfitnesses[maxfit_pop][maxfit_gen,:])[2])  # A random columm with maximum fitness
  prds = preds( mpredecessors[maxfit_pop], maxfit_gen, column )
  [(mphenotypes[maxfit_pop][i,prds[i]],mfitnesses[maxfit_pop][i,prds[i]]) for i = 1:maxfit_gen ]
end

# Returns the number of fitness decreases in the fitness history of a random maximum fitness 
#   element of the first population with maximum fitness.
# Also returns the maximum fitness over all generations and the first generation to reach this fitness.
function preds_fit_decreases( predecessors::Array{Int64,2}, phenotypes::Array{MyInt,2}, fitnesses::Array{Float64,2} )
  #println("preds_fit_decreases: size(phenotypes): ",size(phenotypes))
  ppf = preds_phenos_fits( predecessors, phenotypes, fitnesses )
  maxfit_gen = length(ppf)
  fits = map( x->x[2], ppf )
  decreases = 0
  for i = 1:(maxfit_gen-1)
    if fits[i+1] < fits[i]
      decreases += 1
    end
  end
  (decreases, fits[maxfit_gen], maxfit_gen)
end

function preds_fit_decreases( mpredecessors::Array{Array{Int64,2},1}, mphenotypes::Array{Array{MyInt,2},1}, 
    mfitnesses::Array{Array{Float64,2},1} ) 
  #println("preds_fit_decreases: size(mphenotypes[1]): ",size(mphenotypes[1]))
  ppf = preds_phenos_fits( mpredecessors, mphenotypes, mfitnesses )
  maxfit_gen = length(ppf)
  fits = map( x->x[2], ppf )
  decreases = 0
  for i = 1:(maxfit_gen-1)
    if fits[i+1] < fits[i]
      decreases += 1
    end
  end
  (decreases, fits[maxfit_gen], maxfit_gen)
end

# Return the list of successors of element index of population gen to num_gens
# Usually one will call with gen = 1
function successors( predecessors::Array{Int64,2}, gen::Int64, index::Int64 )
  (num_gens, popsize) = size(predecessors)
  succs = Vector{Int64}[]
  succ = [index]
  for i = gen:num_gens
    succ = findall(x->x in succ,predecessors[i,:])
    push!(succs,succ)
  end
  succs
end

# Check that preds(successors(i) == i (approximately stated)
# Check that for each s in successors(i), preds(s) == i
function test_preds_successors( predecessors::Array{Int64,2} )
  (num_gens, num_indices) = size(predecessors)
  for i = 1:num_indices
    succs = successors( predecessors, 1, i )[num_gens]
    for s in succs
      if i != (prd=preds( predecessors, num_gens, s )[1] )
        println("i: ",i,"  succs: ",succs,"  s: ",s,"  prd: ",prd)
      end
    end
  end
end
    
# Check that for each c=1:popsize, s in successors(preds(c))
function test_successors_preds( predecessors::Array{Int64,2} )
  (num_gens, popsize) = size(predecessors)
  for c = 1:num_indices
    p = preds(predecessors, num_gens, c )[1]
    succs = successors( predecessors, 1, p )[num_gens]
    if !( c in succs )
      println("p:  ",p,"  c: ",c,"  succs: ",succs)
    end
  end
end

# Insert a gate into c so that output_values(c) is unchanged.
# The new gate is inserted at position new_gate_index.
# The following gates are shifted right by 1.
# The new gate replaces an input connection to the existing gate a position new_gate_index (before the shift)
# The output of the new gate is the input to the existing gate.
# One input to the new gate is the source of the input to the exisiting gate. 
# The other is chosen randomly according to the levelsback constraint.
function insert_gate!( c::Chromosome )
  funcs = default_funcs(c.params.numinputs)
  p = c.params
  # Chose a random interior node as the "new gate".
  # Don't insert a gate to replace an output gate:
  new_gate_index = rand(1:(p.numinteriors-p.numoutputs+1))  # Position of the new gate
  println("new_gate_index: ",new_gate_index)
  p = c.params = Parameters( p.numinputs, p.numoutputs, p.numinteriors+1, p.numlevelsback+1 )
  new_gate = deepcopy(c.interiors[new_gate_index])  # will be modified below
  input_index = rand(1:2)
  new_gate.func = input_index == 1 ? IN1 : IN2
  c.interiors = vcat(c.interiors[1:(new_gate_index-1)],[new_gate],c.interiors[new_gate_index:end])
  for i = 1:p.numoutputs
    index = p.numinputs + p.numinteriors + i - p.numoutputs # use the last numoutputs interiors 
    c.outputs[i] = OutputNode(index) 
  end
  for i = (new_gate_index+1):p.numinteriors
    for j = 1:p.nodearity
      if c.interiors[i].inputs[j] >= p.numinputs + new_gate_index 
        c.interiors[i].inputs[j] += 1
      end
    end
  end
  # Modify one of the inputs of interiors[new_gate_index+1] to point to the new gate
  c.interiors[new_gate_index+1].inputs[input_index] = p.numinputs + new_gate_index
  set_active_to_false(c)
  c
end

function delete_gate!( c::Chromosome, interior_to_delete::Int64 )
  #println("delete_gate! interior_to_delete: ",interior_to_delete,"  active: ",c.interiors[interior_to_delete].active)
  p = c.params
  @assert interior_to_delete <= p.numinteriors - p.numoutputs
  gate_to_delete = c.interiors[interior_to_delete]
  new_levelsback = (p.numlevelsback > p.numinputs) ? p.numlevelsback-1 : p.numlevelsback
  p = c.params = Parameters( p.numinputs, p.numoutputs, p.numinteriors-1, new_levelsback )
  c.interiors = vcat(c.interiors[1:(interior_to_delete-1)],c.interiors[(interior_to_delete+1):end])
  for i = 1:p.numoutputs
    index = p.numinputs + p.numinteriors + i - p.numoutputs # use the last numoutputs interiors 
    c.outputs[i] = OutputNode(index) 
  end
  for i = interior_to_delete:p.numinteriors
    for j = 1:p.nodearity
      if c.interiors[i].inputs[j] >= interior_to_delete + p.numinputs
        c.interiors[i].inputs[j] -= 1
      end
    end
  end
  c
end

function delete_gate!( c::Chromosome )
  output_values(c)  # needed to set active gates to active  
  if number_active_gates(c) == c.params.numinteriors
    dg = rand(1:(c.params.numinteriors-p.numoutputs))
    return delete_gate!( c, dg )
  end
  inactive_interior_list = filter!(i->!c.interiors[i].active, collect(1:c.params.numinteriors))
  dg = rand(inactive_interior_list)
  return delete_gate!( c, dg )
end

function test_delete_gate( n::Int64, p::Parameters )
  for i = 1:n
    c = random_chromosome(p)
    sav_c = deepcopy(c)
    output_values(sav_c)  # needed to set active gates to active
    #dg = rand(1:(c.params.numinteriors-p.numoutputs))
    #println("dg: ",dg,"  active: ",sav_c.interiors[dg].active)
    #delete_gate!(c,dg)
    #if sav_c.interiors[dg].active == false
    delete_gate!(c)
    if number_active_gates(sav_c) < sav_c.params.numinteriors
      println("test for i = ",i)
      if output_values(c) != output_values(sav_c)
        println((output_values(c),output_values(sav_c)))
        print_circuit(c)
        print_circuit(sav_c)
        return (sav_c,c)
      end
    end
  end
end 

function test_insert_gate( n::Int64, p::Parameters )
  for i = 1:n
    c = random_chromosome(p)
    sav_c = deepcopy(c)
    insert_gate!(c)
    if output_values(c) != output_values(sav_c)
      println((output_values(c),output_values(sav_c)))
      print_circuit(c)
      print_circuit(sav_c)
      return (sav_c,c)
    end
  end
end 
