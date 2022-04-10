# Circuit entropy, Adami complexity, and mutual information between the genotype and phenotype fitness distributions.
# A dataframe of results is produced by running simple_pop_evolve()
# Test functions are test_update_chrome_values() and check_marginals()
# See data/4_3_22/  for the output of a test run.

# The parameters that describe one chromosome.
# Assumes that gates are arity 2.
# The three components of an interior node (gate) of a chromosome are the function and the two inputs.
# The parameters for the interior node are integers specifying these three components.
function chrome_params( ch::Chromosome, funcs::Vector{Func} )
  arities = gate_arities( p, funcs )
  ch_params = Vector{Int64}[]
  for i = 1:ch.params.numinteriors
    j = 1  # determine index of func of interior node i
    while j <= length(funcs) && ch.interiors[i].func.func != funcs[j].func
      j += 1
    end
    #println("i: ",i,"  j: ",j)
    push!(ch_params, [ j, ch.interiors[i].inputs[1]-arities[i][2][2], ch.interiors[i].inputs[2]-arities[i][3][2] ])
  end
  ch_params
end

# based on random_chromosome() in Chromosome.jl
# The three components of an interior node of a chromosome are the function and the two inputs.
# Each of these is specified by an integer.  
# The corresponding arity is the pair of the number of possible values for the integer and the minimum value of the integer.
# Thus, the arity of an interior node is the triple of these pairs.
# gate_arities() returns the list of arities for each of the interior nodes of a chromosome corresponding to paramters p.
function gate_arities( p::Parameters, funcs::Vector{Func} )
  gate_arities = Vector{Tuple{Int64,Int64}}[]
  for index = (p.numinputs +1):(p.numinputs+p.numinteriors)
    maxindex = index - 1
    minindex = max(1,index-p.numlevelsback) 
    arities = [(length(funcs),1),(maxindex-minindex+1,minindex-1),(maxindex-minindex+1,minindex-1)]
    push!(gate_arities,arities)
  end
  gate_arities
end

# A record of the number of chromosomes that have a value for each of these fields
mutable struct GateValues
  func::Vector{Int64}
  input1::Vector{Int64}
  input2::Vector{Int64}
end

# Returns a list of zero GateValues, one for each of the p.numinteriors interior nodes.
# Each GateValues vector is a zeros vector whose length is specified by arities[1] for that interior node.
function init_chrome_values( p::Parameters, funcs::Vector{Func} )
  arities = gate_arities( p, funcs )
  chrome_values = GateValues[]
  for i = 1:p.numinteriors
    chv = GateValues( zeros(Int64,arities[i][1][1]), zeros(Int64,arities[i][2][1]), zeros(Int64,arities[i][3][1]))
    push!(chrome_values,chv)
  end
  chrome_values
end

# update chrome_values for one chromosome by incrementing one component of the corresponding chrome_values array by 1.
function update_chrome_values!( chrome_values::Vector{GateValues}, arities::Vector{Vector{Tuple{Int64, Int64}}}, ch::Chromosome, funcs::Vector{Func} )
  @assert ch.params.numinteriors <= length(chrome_values)
  #print_circuit(ch)
  for i = 1:ch.params.numinteriors
    j = 1  # determine index of func of interior node i
    while j <= length(funcs) && ch.interiors[i].func.func != funcs[j].func
      j += 1
    end
    #println("i: ",i,"  j: ",j)
    chrome_values[i].func[j] += 1
    chrome_values[i].input1[ch.interiors[i].inputs[1]-arities[i][2][2]] += 1
    chrome_values[i].input2[ch.interiors[i].inputs[2]-arities[i][3][2]] += 1
  end
end

# update chrome_values for a population of chromosomes
function update_chrome_values!( chrome_values::Vector{GateValues}, arities::Vector{Vector{Tuple{Int64, Int64}}}, pop::Vector{Chromosome}, funcs::Vector{Func} )
  for ch in pop
    update_chrome_values!( chrome_values, arities, ch, funcs )
  end
end

# The entropy of a vector of Integers where v[i] is the number of chromosomes with value i for an interior node field
# This is just the standard entropy formula
function aentropy( v::Vector{Int64} )
  log2z(x) = (x == 0.0 ? 0.0 : log2(x))  # log2 with 0 mapped to 0
  sumv = sum(v)
  if sumv != 0 
    p = v/sumv
  else
    return 0.0
  end
  sum( -p[i]*log2z(p[i]) for i = 1:length(p) )
end

function entropy( gv::Vector{GateValues} )
  entropies = Float64[]
  for g in gv
    ents = [aentropy(g.func),aentropy(g.input1),aentropy(g.input2)]
    entropies = vcat(entropies,ents)
  end
  entropies
end
  
function pop_entropies( p::Parameters, funcs::Vector{Func}, pop::Vector{Chromosome} )
  for ch in pop
    @assert ch.params.numinputs == p.numinputs
    @assert ch.params.numoutputs == p.numoutputs
    @assert ch.params.numinteriors == p.numinteriors
    @assert ch.params.numlevelsback == p.numlevelsback
  end
  arities=gate_arities(p,funcs)
  chvalues = init_chrome_values( p, funcs )
  update_chrome_values!( chvalues, arities, pop, funcs )
  #println("chvalues: ",chvalues)
  entropy(chvalues)
end

function avg_pop_entropy( pop::Vector{Chromosome}, funcs::Vector{Func}=default_funcs(pop[1].params) )
  arities=gate_arities(p,funcs)
  chvalues = init_chrome_values( p, funcs )
  update_chrome_values!( chvalues, arities, pop, funcs )
  #println("chvalues: ",chvalues)
  entvalues = entropy(chvalues)
  sum(entvalues)/length(entvalues)
end

function hamming_fitness_nonzero( x::MyInt, target::MyInt, numinputs::Int64 )
  1.0 - hamming_distance( x, target, numinputs)  + 1/2^numinputs
end

# by Standard formula Equation 2.28 of Cover & Thomas
function mut_info( pop::Vector{Chromosome}, target::MyInt, funcs::Vector{Func} ) 
  (pheno_fitness_dict,geno_count_dict,joint_dict) = mut_info_prob_dists( pop, target, funcs )
  sum( joint_dict[kj]*log2(joint_dict[kj]/geno_count_dict[kj[1]]/pheno_fitness_dict[kj[2]]) for kj in keys(joint_dict) )
end

# Creates dictionaries representing fitness probability distributions over phenotypes, genotypes, 
#   and the joint distribution over phenotypes and genoypes.  
# Genotypes that do not occur in the population pop are not included
function mut_info_prob_dists( pop::Vector{Chromosome}, target::MyInt, funcs::Vector{Func} )
  p = pop[1].params
  arities=gate_arities(p,funcs)    
  pheno_fitness_dict = Dict{MyInt,Float64}()
  geno_count_dict = Dict{Vector{Vector{Int64}},Float64}()
  joint_dict = Dict{Tuple{Vector{Vector{Int64}},MyInt},Float64}()
  for ch in pop
    pheno = output_values( ch )[1]
    ph_fit = hamming_fitness_nonzero( pheno, target, p.numinputs )
    pheno_fitness_dict[pheno] = get( pheno_fitness_dict, pheno, 0.0 ) + ph_fit 
    ch_params = chrome_params( ch, funcs )
    geno_count_dict[ch_params] = get( geno_count_dict, ch_params, 0.0 ) + ph_fit
    joint_dict[(ch_params,pheno)] = get( joint_dict, (ch_params,pheno), 0.0 ) + ph_fit
  end
  ph_sum = sum( pheno_fitness_dict[k] for k in keys(pheno_fitness_dict) ) 
  for kp in keys(pheno_fitness_dict)
    pheno_fitness_dict[kp] = pheno_fitness_dict[kp]/ph_sum
  end
  geno_sum = sum( geno_count_dict[k] for k in keys(geno_count_dict) ) 
  for kg in keys(geno_count_dict)
    geno_count_dict[kg] = geno_count_dict[kg]/geno_sum
  end
  #(pheno_fitness_dict, geno_count_dict)
  joint_sum = sum( joint_dict[k] for k in keys(joint_dict) )
  for kj in keys(joint_dict)
    joint_dict[kj] = joint_dict[kj]/joint_sum
  end
  (pheno_fitness_dict,geno_count_dict,joint_dict)
end

# Computes marginal vectors for the genotype and phenotype distributions computed by mut_info_prob_dists()
function marginals( joint_dict::Dict{Tuple{Vector{Vector{Int64}},MyInt},Float64}, 
    pheno_fitness_dict::Dict{MyInt,Float64}, geno_count_dict::Dict{Vector{Vector{Int64}},Float64} )
  geno_marginals = [ sum( get( joint_dict, (kg,kf), 0.0 ) for kf in keys(pheno_fitness_dict) ) for kg in keys(geno_count_dict) ]
  pheno_marginals = [ sum( get( joint_dict, (kg,kf), 0.0 ) for kg in keys(geno_count_dict) ) for kf in keys(pheno_fitness_dict) ] 
  (geno_marginals,pheno_marginals)
end

# Checks that the marginal vectors agree with the genotype and phenotype distributions computed by mut_info_prob_dists() 
# An important check of correctness.
function check_marginals( joint_dict::Dict{Tuple{Vector{Vector{Int64}},MyInt},Float64},
    pheno_fitness_dict::Dict{MyInt,Float64}, geno_count_dict::Dict{Vector{Vector{Int64}},Float64} )
  (gm,pm) = marginals( joint_dict, pheno_fitness_dict, geno_count_dict )
  gdist = [geno_count_dict[k] for k in keys(geno_count_dict)]
  pdist = [pheno_fitness_dict[k] for k in keys(pheno_fitness_dict)]
  @assert gm==gdist
  @assert pm==pdist
  [ (gm,gdist),(pm,pdist) ]
end
  
# Tests only that there are no indexing errors
function test_update_chrome_values(p,popsize::Int64,funcs=default_funcs(p))
  pop = Chromosome[]
  for i = 1:popsize
    push!(pop,random_chromosome(p,funcs))
  end
  arities=gate_arities(p,funcs)
  chvalues = init_chrome_values( p, funcs )
  update_chrome_values!( chvalues, arities, pop, funcs )
  chvalues
end

#=
p = Parameters(2,1,4,3)
funcs = default_funcs(p)
println("arities: ",circuit_arities(p,funcs))
popsize = 8
pop = [ random_chromosome(p,funcs) for _=1:popsize ];
(ph,gn,jt) = mut_info(pop,target,funcs); 
check_marginals( jt, ph, gn );
circuit_pop_entropy( p, funcs, pop )
for j = 1:popsize
  println(ch_values(pop[j],3 ))
end
=#
# A simple no-frills pop_evolve().
# Copied from simple_pop_evolve.jl into adami.jl on 4/3/22

function random_population( p::Parameters, popsize::Int64, funcs::Vector{Func} )
  pop = [ random_chromosome(p,funcs) for i = 1: popsize ]
end

# Does a population genetic algorithm and creates a simple dataframe including avg_pop_entropy
#   and mutual_information between genotype and phenotype distributions
function simple_pop_evolve( pop::Vector{Chromosome}, gl::GoalList, ngens::Int64, mutrate::Float64=1.0;
    prdebug::Bool=false )
  df = DataFrame()
  df.gen = Int64[]
  df.fract_optimal = Float64[]
  df.avg_pop_entropy = Float64[]
  df.mutual_inf = Float64[]
  funcs = default_funcs( p.numinputs )
  target = gl[1][1]
  for gen = 1:ngens
    fitness_vector = rescale_fitnesses([ pop[i].fitness for i = 1:popsize ])
    prdebug ? println("fit_vect: ",fitness_vector) : nothing
    prdebug ? println("gen: ",gen,"  max fit: ",maximum(fitness_vector)) : nothing
    propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] )
    prdebug ? println("after propsel gen: ",gen) : nothing
    for c in pop
      #c.fitness = hamming_distance( gl, output_values(c), p.numinputs )
      c.fitness = fitness_funct( p, c, gl )
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=false,include_pheno=true) : nothing
    end       
    #println("fract opt: ",@sprintf("%3.2f",fract_optimal_chromes( pop )),"  ave_ent: ",@sprintf("%3.2f",avg_pop_entropy( pop )))
    push!(df,[gen,fract_optimal_chromes(pop),avg_pop_entropy(pop),mut_info(pop,target,funcs)])
    sav_pop = deepcopy( pop )  
    for i in 1:popsize
      c = pop[i]
      if rand() <= mutrate
        pop[i] = c = mutate_chromosome!(deepcopy(c),funcs)[1]   
      end
    end
  end
  #pop
  df
end

function print_pop( pop::Vector{Chromosome} )
  for c in pop
    print_circuit(c,include_fitness=true,include_robustness=false,include_pheno=true)
  end
end

function fract_optimal_chromes( pop::Vector{Chromosome} )
  count = 0
  for c in pop
    if c.fitness == 1.0
      count += 1
    end
  end
  count/length(pop)
end

# Copied from Pop_evolve.jl.  Remove
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

# Copied from Pop_evolve.jl.  Remove
function fitness_funct( p::Parameters, c::Chromosome, gl::GoalList )
  maximum( 1.0-hamming_distance( g, output_values(c), p.numinputs ) for g in gl )
end
