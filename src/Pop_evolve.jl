# Population-based evolution to try to approximately replicate Labar & Adami 2016.
# The Julia type PredType is defined in Chromosome.jl
using Main.CGP
using Statistics

function run_both_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, goallist::GoalList, ngens::Int64, numpops::Int64;
    csvfile::String="", uniform_start::Bool=true, prdebug::Bool=false )
  df=DataFrame()
  df.goal = Goal[]
  df.ngens = Int64[]
  df.num_pops = Int64[]
  df.popsize = Int64[]
  df.evolvability = Int64[]
  df.maxfit = Float64[]
  df.maxfit_gen = Int64[]
  df_list = pmap(goal->run_both_one_goal(nreps,p,popsize,goal,ngens,numpops,uniform_start=uniform_start),goallist)
  df = vcat( df_list... )
  #=
  for goal in goallist
    println("goal: ",goal)
    df = run_pop_evolve(nreps,p,popsize,goal,ngens,df=df,uniform_start=uniform_start)
    mpopsize = Int64(round(popsize/numpops))   # popsize should be divisible by numpops
    df = run_multiple_pop_evolve(nreps,p,mpopsize,goal,ngens,numpops,df=df,uniform_start=uniform_start)
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

function run_both_one_goal( nreps::Int64, p::Parameters, popsize::Int64, goal::Goal, ngens::Int64, numpops::Int64;
    uniform_start::Bool=true, prdebug::Bool=false )
  df = DataFrame()
  df.goal = Goal[]
  df.ngens = Int64[]
  df.num_pops = Int64[]
  df.popsize = Int64[]
  df.evolvability = Int64[]
  df.maxfit = Float64[]
  df.maxfit_gen = Int64[]
  println("goal: ",goal)
  df = run_pop_evolve(nreps,p,popsize,goal,ngens,df=df,uniform_start=uniform_start)
  mpopsize = Int64(round(popsize/numpops))   # popsize should be divisible by numpops
  df = run_multiple_pop_evolve(nreps,p,mpopsize,goal,ngens,numpops,df=df,uniform_start=uniform_start)  
end

function run_multiple_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, g::Goal, ngens::Int64, numpops::Int64;
    df::DataFrame=DataFrame(), csvfile::String="", uniform_start::Bool=true, prdebug::Bool=false )
  if size(df)[1] == 0
    df.goal = Goal[]
    df.ngens = Int64[]
    df.num_pops = Int64[]
    df.popsize = Int64[]
    df.evolvability = Int64[]
    df.maxfit = Float64[]
    df.maxfit_gen = Int64[]
  end
  for rep = 1:nreps
    (mphenotypes,mpredecessors,mfitnesses,mcomplexities,mnumgates)=multiple_pop_evolve(p,popsize,g,ngens,numpops,uniform_start=false)
    push!(df.goal,g)
    push!(df.ngens,ngens)
    push!(df.num_pops,numpops)
    push!(df.popsize,popsize)
    push!(df.evolvability,final_evolvability(mphenotypes))
    findmax_mfh=findmax(max_fitness_history(mfitnesses))
    push!(df.maxfit,findmax_mfh[1])
    push!(df.maxfit_gen,findmax_mfh[2])
  end
  df
end
    
function run_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, g::Goal, ngens::Int64;
    df::DataFrame=DataFrame(), csvfile::String="", uniform_start::Bool=true, prdebug::Bool=false )
  if size(df)[1] == 0
    df.goal = Goal[]
    df.ngens = Int64[]
    df.num_pops = Int64[]
    df.popsize = Int64[]
    df.evolvability = Int64[]
    df.maxfit = Float64[]
    df.maxfit_gen = Int64[]
  end
  for rep = 1:nreps
    (pop,phenotypes,predecessors,fitnesses,complexities,numgates)=pop_evolve(p,popsize,g,ngens,uniform_start=false)
    push!(df.goal,g)
    push!(df.ngens,ngens)
    push!(df.num_pops,1)
    push!(df.popsize,popsize)
    push!(df.evolvability,final_evolvability(phenotypes))
    findmax_mfh=findmax(max_fitness_history(fitnesses))
    push!(df.maxfit,findmax_mfh[1])
    push!(df.maxfit_gen,findmax_mfh[2])
  end
  df
end
    
function multiple_pop_evolve( p::Parameters, popsize::Int64, g::Goal, ngens::Int64, numpops::Int64;
    uniform_start::Bool=true, prdebug::Bool=false )
  mphenotypes = fill(zeros(MyInt,0,0),numpops)
  mpredecessors = fill(zeros(PredType,0,0),numpops)
  mfitnesses = fill(zeros(Float64,0,0),numpops)  
  mcomplexities = fill(zeros(Float64,0,0),numpops)
  mnumgates = fill(zeros(Int64,0,0),numpops)
  for pind = 1:numpops
    #println("pind: ",pind)
    (pop,mphenotypes[pind],mpredecessors[pind],mfitnesses[pind],mcomplexities[pind],mnumgates[pind]) = 
        pop_evolve(p,popsize,g,ngens,uniform_start=false)
  end
  (mphenotypes,mpredecessors,mfitnesses,mcomplexities,mnumgates)
end

function pop_evolve( p::Parameters, popsize::Int64, g::Goal, maxgenerations::Int64; 
    uniform_start::Bool=true, prdebug::Bool=false )
  funcs = default_funcs( p.numinputs )
  phenotypes = zeros(MyInt,maxgenerations,popsize)
  predecessors = zeros(PredType,maxgenerations,popsize)
  fitnesses = zeros(Float64,maxgenerations,popsize)
  complexities = zeros(Float64,maxgenerations,popsize)
  numgates = zeros(Int64,maxgenerations,popsize)
  if uniform_start
    c = random_chromosome(p,funcs,ident=PredType(1))
    pop = fill( c, popsize )
  else
    pop = [ random_chromosome(p,funcs,ident=PredType(i)) for i = 1: popsize ]
  end
  for i in 1:popsize
    c = pop[i]
    c.robustness = i
    c.fitness = hamming_distance( g, output_values(c), p.numinputs )
    prdebug ? print(c.robustness,"  ") : nothing
    prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
  end
  for gen = 1:maxgenerations
    fitness_vector = rescale_fitnesses([ pop[i].fitness for i = 1:popsize ])
    prdebug ? println("fit_vect: ",fitness_vector) : nothing
    prdebug ? println("gen: ",gen,"  max fit: ",maximum(fitness_vector)) : nothing
    propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] )
    prdebug ? println("after propsel gen: ",gen) : nothing
    for c in pop
      c.fitness = hamming_distance( g, output_values(c), p.numinputs )
      prdebug ? print(c.robustness,"  ") : nothing
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
    end
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
    sav_pop = deepcopy( pop )
    pop = map(c->mutate_chromosome!(deepcopy(c),funcs)[1], pop )
    for i in 1:popsize
      c = pop[i]
      c.robustness = i  
      c.fitness = hamming_distance( g, output_values(c), p.numinputs )
    end
    prdebug ? println("after mutate gen: ",gen) : nothing
    for i in 1:popsize
      c = pop[i]
      c.robustness = i
      c.fitness = hamming_distance( g, output_values(c), p.numinputs )
      prdebug ? print(c.robustness,"  ") : nothing
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
    end
  end
  (pop,phenotypes,predecessors,fitnesses,complexities,numgates)
end

function rescale_fitnesses( fit_vect::Vector{Float64} )
  #fit_min = minimum( fit_vect )
  fit_min = quantile( fit_vect, 0.20 )
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
  history = [maximum( histories[i][g] for i = 1:numpops ) for g = 1:num_gens ]
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
