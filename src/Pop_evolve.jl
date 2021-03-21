# Population-based evolution to try to approximately replicate Labar & Adami 2016.
# The Julia type PredType is defined in Chromosome.jl
using Main.CGP
using Statistics
using DataFramesMeta
using HypothesisTests

# Run both one population and numpops subpopulations to eithe return a dataframe or
#   write the dataframe to csvfile if it is given as a keyword argument
# mutrate is the probability that a chromosome (circuit) will be mutation in a generation
function run_both_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, goallistlist::Vector{GoalList}, 
    ngens::Int64, numpops::Int64, mutrate::CGP.FloatRange=1.0;
    csvfile::String="", uniform_start::Bool=false, prdebug::Bool=false )
  df=DataFrame()
  df.goallist = GoalList[]
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numgates = Int64[]
  df.levelsback = Int64[]
  df.ngens = Int64[]
  df.numpops = Int64[]
  df.popsize = Int64[]
  df.mutrate = Float64[]
  df.evolvability = Int64[]
  df.maxfit = Float64[]
  df.maxfit_gen = Int64[]
  df.fit_decreases = Int64[]
  
  mr_goallists_list = [ (mutr,gl) for mutr = mutrate for gl in goallistlist ]
  df_list = pmap(mr_goallist_pair->run_both_one_goallist(nreps,p,popsize,ngens,numpops,mr_goallist_pair,uniform_start=uniform_start),mr_goallists_list)
  #df_list = map(mr_goallist_pair->run_both_one_goallist(nreps,p,popsize,ngens,numpops,mr_goallist_pair,uniform_start=uniform_start),mr_goallists_list)
  df = vcat( df_list... )
  
  #=
  for goallist in goallistlist
    println("goallist: ",goallist)
    df = run_pop_evolve(nreps,p,popsize,goallist,ngens,df=df,uniform_start=uniform_start)
    mpopsize = Int64(round(popsize/numpops))   # popsize should be divisible by numpops
    df = run_multiple_pop_evolve(nreps,p,mpopsize,goallist,ngens,numpops,mutrate,df=df,uniform_start=uniform_start)
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
      println(f,"# uniform start: ",uniform_start)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Helper function for run_both_pop_evolve().  Used in pmap()
function run_both_one_goallist( nreps::Int64, p::Parameters, popsize::Int64, ngens::Int64, numpops::Int64,
    mr_goal_pair; uniform_start::Bool=false, prdebug::Bool=false )
  (mutrate,goal) = mr_goal_pair 
  df = DataFrame()
  df.goallist = GoalList[]
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numgates = Int64[]
  df.levelsback = Int64[]
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

function run_multiple_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, numpops::Int64,
    mutrate::Float64=1.0;
    df::DataFrame=DataFrame(), csvfile::String="", uniform_start::Bool=false, prdebug::Bool=false )
  println("run_multiple_pop_evolve  mutrate: ",mutrate,"   goallist: ",gl)
  if size(df)[1] == 0
    df.goallist = GoalList[]
    df.numinputs = Int64[]
    df.numoutputs = Int64[]
    df.numgates = Int64[]
    df.levelsback = Int64[]
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
    (mutrate,mphenotypes,mpredecessors,mfitnesses,mcomplexities,mnumgates)=multiple_pop_evolve(p,popsize,gl,ngens,numpops,mutrate,uniform_start=uniform_start)
    #println("run_multiple_pop_evolve size(mphenotypes): ",size(mphenotypes),"  size(mphenotypes[1]): ",size(mphenotypes[1]))
    push!(df.goallist,gl)
    push!(df.numinputs,p.numinputs)
    push!(df.numoutputs,p.numoutputs)
    push!(df.numgates,p.numinteriors)
    push!(df.levelsback,p.numlevelsback)
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
 
function run_pop_evolve( nreps::Int64, p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, mutrate::Float64=1.0;
    df::DataFrame=DataFrame(), csvfile::String="", uniform_start::Bool=false, prdebug::Bool=false )
  println("run_pop_evolve  mutrate: ",mutrate,"   goallist: ",gl)
  if size(df)[1] == 0
    df.goallist = GoalList[]
    df.numinputs = Int64[]
    df.numoutputs = Int64[]
    df.numgates = Int64[]
    df.levelsback = Int64[]
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
        pop_evolve(p,popsize,gl,ngens,mutrate,uniform_start=uniform_start)
    #println("run_pop_evolve size(phenotypes): ",size(phenotypes))
    push!(df.goallist,gl)
    push!(df.numinputs,p.numinputs)
    push!(df.numoutputs,p.numoutputs)
    push!(df.numgates,p.numinteriors)
    push!(df.levelsback,p.numlevelsback)
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
    
function multiple_pop_evolve( p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, numpops::Int64, mutrate::Float64=1.0;
    uniform_start::Bool=false, prdebug::Bool=false )
  mphenotypes = fill(zeros(MyInt,0,0),numpops)
  mpredecessors = fill(zeros(PredType,0,0),numpops)
  mfitnesses = fill(zeros(Float64,0,0),numpops)  
  mcomplexities = fill(zeros(Float64,0,0),numpops)
  mnumgates = fill(zeros(Int64,0,0),numpops)
  for pind = 1:numpops
    #println("pind: ",pind)
    (pop,mphenotypes[pind],mpredecessors[pind],mfitnesses[pind],mcomplexities[pind],mnumgates[pind]) = 
        pop_evolve(p,popsize,gl,ngens,mutrate,uniform_start=false)
  end
  (mutrate,mphenotypes,mpredecessors,mfitnesses,mcomplexities,mnumgates)
end

function pop_evolve( p::Parameters, popsize::Int64, gl::GoalList, ngens::Int64, mutrate::Float64=1.0; 
    uniform_start::Bool=false, prdebug::Bool=false )
  #println("pop_evolve: goals: ",gl,"  psize: ",popsize,"  ngens: ",ngens,"  mutrate: ",mutrate)
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
    c.fitness = fitness_funct( p, c, gl )
    prdebug ? print(c.robustness,"  ") : nothing
    prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
  end
  for gen = 1:ngens
    fitness_vector = rescale_fitnesses([ pop[i].fitness for i = 1:popsize ])
    #fitness_vector = ([ pop[i].fitness for i = 1:popsize ])
    prdebug ? println("fit_vect: ",fitness_vector) : nothing
    prdebug ? println("gen: ",gen,"  max fit: ",maximum(fitness_vector)) : nothing
    propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] )
    prdebug ? println("after propsel gen: ",gen) : nothing
    for c in pop
      #c.fitness = hamming_distance( gl, output_values(c), p.numinputs )
      c.fitness = fitness_funct( p, c, gl )
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
        #c.fitness = 1.0-hamming_distance( g, output_values(c), p.numinputs )
        c.fitness = fitness_funct( p, c, gl )
        #println("mutation count: ",count_mutations,"  phenotype: ",@sprintf("0x%x",(output_values(c)[1])),"  fitness: ",c.fitness)
        count_mutations += 1
      end
      c.robustness = i  # Robustness is not robustness; instead it is an identifier of the index of circut in the population
    end
    
    for i in 1:length(pop)
      phenotypes[gen,i] = output_values(pop[i])[1]
      predecessors[gen,i] = pop[i].robustness
      fitnesses[gen,i] = pop[i].fitness = fitness_funct( p, c, gl ) 
      @assert pop[i].fitness == fitness_funct( p, c, gl )
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
      #c.fitness = 1.0-hamming_distance( g, output_values(c), p.numinputs )
      c.fitness = fitness_funct( p, c, gl )
      prdebug ? print(c.robustness,"  ") : nothing
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=true,include_pheno=true) : nothing
    end
    =#
  end
  #println("end pop evolve: len(pop): ",length(pop),"  size: ",size(phenotypes))
  (pop,phenotypes,predecessors,fitnesses,complexities,numgates)
end

function fitness_funct( p::Parameters, c::Chromosome, gl::GoalList )
  maximum( 1.0-hamming_distance( g, output_values(c), p.numinputs ) for g in gl )
  #complexity5(c)
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

function randgoallistlist( ngoals::Int64, ngoallists::Int64, p::Parameters )
  [ randgoallist( ngoals, p ) for _ = 1:ngoallists ]
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

function compare_fract_decreases( df::DataFrame, csvfile::String )
  cdf=compare_fract_decreases( df )
  write_dataframe_with_comments(cdf,csvfile,"$(csvfile[1:(end-4)])fract_dec.csv")
  cdf 
end

function compare_fract_decreases( csvfile::String )
  df = read_dataframe( csvfile )
  cdf=compare_fract_decreases( df )
  write_dataframe_with_comments(cdf,csvfile,"$(csvfile[1:(end-4)])fract_dec.csv")
  cdf 
end

# Returns a triple of the large population mean, the small population mean, and the p-value of the Mann Whitney test
function compare_fract_decreases( df::DataFrame )
  df.fract_dec = df.fit_decreases./df.maxfit_gen
  num_small_pops = unique(df.numpops)[2]    # The number of pops in the small pop case
  cdf=DataFrame()
  if ("numinputs" in names(df)) || (:numinputs in names(df))
    cdf.numinputs = Int64[]
    cdf.numoutputs = Int64[]
    cdf.numgates = Int64[]
    cdf.levelsback = Int64[]
  end
  cdf.ngens=Int64[]
  cdf.lg_popsize = Int64[]
  cdf.sm_popsize = Int64[]
  cdf.mutrate = Float64[]
  cdf.lge_pop_fdecr = Float64[]
  cdf.sml_pop_fdecr = Float64[]
  cdf.MannWhitPval = Float64[]
  for mr in unique(df.mutrate)
    ldf = @where(df,(:maxfit.>=1.0).&(:numpops.==1).&(:mutrate.==mr))
    large_pop_fract_dec = mean(ldf.fract_dec)   
    sdf = @where(df,(:maxfit.>=1.0).&(:numpops.==num_small_pops).&(:mutrate.==mr))
    small_pops_fract_dec = mean(sdf.fract_dec)   
    println("mr: ",mr,"  length(ldf.fract_dec): ",length(ldf.fract_dec),"  length(sdf.fract_dec): ",length(sdf.fract_dec))
    if length( sdf.fract_dec ) == 0 || length( ldf.fract_dec ) == 0
      println("warning: length(sdf) == 0  or length(ldf) == 0 for mutrate: ",mr)
    else
      MWT = pvalue(MannWhitneyUTest( ldf.fract_dec, sdf.fract_dec )) 
      if ("numinputs" in names(df)) || (:numinputs in names(df))
        push!(cdf,(sdf.numinputs[1],sdf.numoutputs[1],sdf.numgates[1],sdf.levelsback[1],
            sdf.ngens[1],ldf.popsize[1],sdf.popsize[1],mr,large_pop_fract_dec,small_pops_fract_dec,MWT))
      else
        push!(cdf,(sdf.ngens[1],ldf.popsize[1],sdf.popsize[1],mr,large_pop_fract_dec,small_pops_fract_dec,MWT))
      end
    end
  end
  cdf
end
