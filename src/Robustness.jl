export robustness, genotype_robustness, phenotype_robustness, matrix_robustness

# Computes the genotype robustness of circuit c.
function robustness( c::Circuit, funcs::Vector{Func} )
  #print("robustness: c:  ")
  #print_circuit(c,funcs)
  c_output = output_values(c,funcs)
  outputs = mutate_all( c, funcs, output_outputs=true )
  #println("outputs[1]: ",outputs[1])
  robust_outputs = filter( x->x==c_output, outputs )
  return length(robust_outputs)/length(outputs)
end   

function grobustness( c::Circuit, funcs::Vector{Func} )
  #print("robustness: c:  ")
  #print_circuit(c,funcs)
  c_output = output_values(c,funcs)
  outputs = mutate_all( c, funcs, output_outputs=true )
  #println("outputs[1]: ",outputs[1])
  robust_outputs = filter( x->x==c_output, outputs )
  return length(robust_outputs)/length(outputs)
end   

# Returns a list of the robustnesses of numcircuits random genotypes
function genotype_robustness( p::Parameters, funcs::Vector{Func}, numcircuits::Int64; use_lincircuit::Bool=false )
  robust_list = Float64[]
  for i = 1:numcircuits
    c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
    push!( robust_list, robustness( c, funcs ) )
  end
  robust_list
end

# If summarize==false
#   returns a list of the (robustness,evolvability) pairs for numcircuits evolved genotypes with the given parameters and funcs
# If summarize==true
#   returns the (mean(robustness),mean(evolvability0) pair for numcircuits evolved genotypes with the given parameters and funcs
function genotype_robustness_evol( p::Parameters, funcs::Vector{Func}, pheno::Goal, numcircuits::Int64, max_steps::Int64=200_000, max_tries::Int64=10; 
      use_lincircuit::Bool=false, summarize::Bool=false )
  robust_evo_list = Tuple{Float64,Float64}[]
  for i = 1:numcircuits
    #c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
    c = pheno_evolve( p, funcs, pheno, max_tries, max_steps )[1]
    if c != nothing
      push!( robust_evo_list, ( robustness( c, funcs ), genotype_evolvability( c, funcs ) ) )
    else
      continue
    end
  end
  if summarize 
    (mean(map(x->x[1],robust_evo_list)),mean(map(x->x[2],robust_evo_list)))
  else
    robust_evo_list
  end
end

function genotype_robustness_evol( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, numcircuits::Int64, max_steps::Int64=200_000, max_tries::Int64=10; 
      use_lincircuit::Bool=false, summarize::Bool=false )
  println("function genotype_robustness_evol(): inputs: ",p.numinputs,"  gates: ",p.numinteriors )
  df = DataFrame( :goal=>Goal[], :robust=>Float64[], :evolvability=>Float64[] )
  #robust_list = Float64[]
  #evolvability_list = Float64[]
  for ph in ph_list
    ( robust, evolvability ) = genotype_robustness_evol( p14, funcs, ph, 8, summarize=true )
    push!( df, ( ph, robust, evolvability ) )
  end
  df
end

# Average robustness of numcircuits random chromosomes for a list of parameters settings
# Triple:  (numiinputs, numinteriors, numlevelsback)  
# Example:  genotype_robustness( [(3,8,4),(4,10,5)], 100 )
# Note there is another function with this name above.
function genotype_robustness( ni_ng_lb_triples::Vector{Tuple{Int64,Int64,Int64}}, numcircuits::Int64 )
  ni_list = Int64[]
  ng_list = Int64[]
  lb_list = Int64[]
  mean_rbst_list = Float64[]
  std_rbst_list = Float64[]
  q10_rbst_list = Float64[]
  q90_rbst_list = Float64[]
  rlist_list = Vector{Float64}[]
  for (ni, ng, lb) in ni_ng_lb_triples
    p = Parameters( ni, 1, ng, lb )
    funcs = default_funcs( p )
    rlist = Float64[]
    for i = 1:numcircuits
      push!(rlist,robustness(random_chromosome(p,funcs),funcs))
    end
    push!( ni_list, ni )
    push!( ng_list, ng )
    push!( lb_list, lb )
    push!( mean_rbst_list, mean(rlist ) )
    push!( std_rbst_list, std(rlist ) )
    push!( q10_rbst_list, quantile( rlist, 0.1 ) )
    push!( q90_rbst_list, quantile( rlist, 0.9 ) )
    push!( rlist_list, rlist )
  end
  df = DataFrame( :ninteriors=>ni_list, :ngates=>ng_list, :lb=>lb_list, :numcircuits=>fill(nreps,length(ni_ng_lb_triples)), 
      :mean_robust=>mean_rbst_list, :std_robust=>std_rbst_list, :q10_robust=>q10_rbst_list, :q90_robust=>q90_rbst_list,
      :rlists=>rlist_list )
end

# Alternative: function run_ph_evolve() in Evolve.jl
function phenotype_robustness( ni_ng_lb_triples::Vector{Tuple{Int64,Int64,Int64}}, numcircuits::Int64, ngoals::Int64, max_tries::Int64, max_steps::Int64;
     use_lincircuit::Bool=false, use_mut_evolve::Bool=false, csvfile::String="" )
  ni_list = Int64[]
  ng_list = Int64[]
  lb_list = Int64[]
  mean_rbst_list = Float64[]
  std_rbst_list = Float64[]
  q10_rbst_list = Float64[]
  q90_rbst_list = Float64[]
  rlist_list = Vector{Float64}[]
  funcs = nothing  # establish scope
  for (ni, ng, lb) in ni_ng_lb_triples
    p = Parameters( ni, 1, ng, lb )
    funcs = default_funcs( p )
    println("funcs: ",funcs)
    phlist = sort(randgoallist(ngoals,p))
    rdf = run_pheno_evolve_rbst( p, funcs, phlist, numcircuits, max_tries, max_steps )
    push!( ni_list, ni )
    push!( ng_list, ng )
    push!( lb_list, lb )
    push!( mean_rbst_list, mean( rdf.mean_rbst) )
    push!( std_rbst_list, mean( rdf.std_rbst ) )
    push!( q10_rbst_list, mean( rdf.q10_rbst ) )
    push!( q90_rbst_list, mean( rdf.q90_rbst ) )
  end
  df = DataFrame( :numinputs=>ni_list, :ngates=>ng_list, :lb=>lb_list, :nfuncs=>length(funcs), :numcircuits=>fill(numcircuits,length(ni_ng_lb_triples)), 
        :ngoals=>fill(ngoals,length(ni_ng_lb_triples)), :mean_robust=>mean_rbst_list, :std_robust=>std_rbst_list, :q10_robust=>q10_rbst_list, :q90_robust=>q90_rbst_list )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# use_lincircuit: ",use_lincircuit)
      println(f,"# MyInt: ",MyInt)
      println(f,"# ni_ng_lb_triples: ",ni_ng_lb_triples)
      println(f,"# ngoals: ",ngoals)
      println(f,"# numcircuits: ",numcircuits)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write(f, df, append=true, writeheader=true )
    end
  end
  df
end

# Robustness for 1 phenotype
function matrix_robustness( p::Parameters, funcs::Vector{Func}, E::Matrix, ph::Goal )
  @assert 2^2^p.numinputs == size(E)[1]
  E[ph[1]+1,ph[1]+1]/sum(E[ph[1]+1,:])
end

# Robustness for a list of phenotypes
function matrix_robustness( p::Parameters, funcs::Vector{Func}, E::Matrix, phlist::GoalList )
  @assert 2^2^p.numinputs == size(E)[1]
  map( ph->matrix_robustness( p, funcs, E, ph ), phlist )
end

# Evolves numciruits circuits to map to each phenotype in phlist.
# Returns a dataframe containing robustness statistics and rlists.
# Alternative: function run_ph_evolve() in Evolve.jl
function run_pheno_evolve_rbst( p::Parameters, funcs::Vector{Func}, phlist::GoalList, numcircuits::Int64, max_tries::Int64, max_steps::Int64;
    use_lincircuit::Bool=false, use_mut_evolve::Bool=false, csvfile::String="" )
  nphenos = length(phlist)
  if max_tries <= numcircuits
    error("function run_ph_evolve(): max_tries should be greater than numcircuits. Your values: max_tries: ",max_tries,"  numcircuits: ",numcircuits)
  end
  #rdict = redundancy_dict(p,funcs)
  #kdict = kolmogorov_complexity_dict(p,funcs)
  nphenos = length(phlist)
  # result_list is a list of (circuit,steps) pairs
  #result_list = map( ph->pheno_evolve( p, funcs, ph, numcircuits, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve), phlist ) 
  result_list = pmap( ph->pheno_evolve( p, funcs, ph, numcircuits, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve), phlist ) 
  mean_rbst_list = Float64[]
  std_rbst_list = Float64[]
  q10_rbst_list = Float64[]
  q90_rbst_list = Float64[]
  rlist_list = Vector{Float64}[]
  #println("result_list: ",result_list)
  for res in result_list
    nc_list = map( r->r[1], res )   # circuits from res
    rlist = Float64[]
    for nc in nc_list
      push!(rlist,robustness(nc,funcs))
    end  
    push!( mean_rbst_list, mean(rlist ) )
    push!( std_rbst_list, std(rlist ) )
    push!( q10_rbst_list, quantile( rlist, 0.1 ) )
    push!( q90_rbst_list, quantile( rlist, 0.9 ) )
    push!( rlist_list, rlist )
  end
  df = DataFrame( :pheno=>phlist, :mean_rbst=>mean_rbst_list, :std_rbst=>std_rbst_list, 
    :q10_rbst=>q10_rbst_list, :q90_rbst=>q90_rbst_list, :rlists=>rlist_list )
end
#=
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# MyInt: ",MyInt)
      println(f,"# ",use_lincircuit ? "LGP" : "CGP" )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# length(phlist): ",length(phlist))
      println(f,"# numcircuits: ",numcircuits)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
=#

# jupyter notebook:  data/2_26_23/pheno_robust2_26_23.ipynb

# Phenotype robustness under the random model of Greenbury "Genetic Correlations" 2016
function random_pheno_robustness( p::Parameters, funcs::Vector{Func}, ph::Goal, rdict::Dict=redundancy_dict(p,funcs) )
  redund = rdict[ph[1]]
end

# Phenotype robustness computed based the circ_ints_list of a counts dataframe.
# Example:  p = Parameters( 3, 1, 8, 4 ); funcs=default_funcs(p)[1:4]
#     wdf = read_dataframe("../data/counts/count_outputs_ch_4funcs_3inputs_8gates_4lb_ge_W.csv")
#     sampling_pheno_robustness( p, funcs, [0x000d], wdf.circuits_list )
function sampling_pheno_robustness( p::Parameters, funcs::Vector{Func}, ph::Goal, circ_ints_list::Union{Vector{Vector{Int128}},Vector{String}} )
  if typeof(circ_ints_list) == Vector{String}
    circ_ints = string_to_expression( circ_ints_list[ph[1]+1] )
  else
    circ_ints = circ_ints_list[ph[1]+1]
  end
  @assert output_values( int_to_chromosome( circ_ints[1], p, funcs )) == ph
  circuits = map( ci->int_to_chromosome( ci, p, funcs ), circ_ints )
  mean( map( ch->robustness( ch, funcs ), circuits ) )
end

# Phenotype robustnesses computed based the circ_ints_list of a counts dataframe.
# Example:  sampling_pheno_robustnesses( p, funcs, [[0x000d],[0x0033]], wdf.circuits_list )
function sampling_pheno_robustnesses( p::Parameters, funcs::Vector{Func}, phlist::GoalList, circ_ints_list::Union{Vector{Vector{Int128}},Vector{String}} )
  map( ph->sampling_pheno_robustness( p, funcs, ph, circ_ints_list), phlist )
end

# Phenotype robustness computed based evolving ncircuits circuits that map to ph
# Examle:  evolution_pheno_robustness( p, funcs, [0x000d], 5, 10, 100_000 )
function evolution_pheno_robustness( p::Parameters, funcs::Vector{Func}, ph::Goal, ncircuits::Int64, max_tries::Int64, max_steps::Int64 )
  circ_steps_list = pheno_evolve( p, funcs, ph, ncircuits, max_tries, max_steps )
  circuits_list = map( x->x[1], circ_steps_list )
  mean( map( ch->robustness( ch, funcs ), circuits_list ) )
end

# Phenotype robustnesses computed based evolving ncircuits circuits that map to each ph in phlist
function evolution_pheno_robustnesses( p::Parameters, funcs::Vector{Func}, phlist::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64 )
  map( ph->evolution_pheno_robustness( p, funcs, ph, ncircuits, max_tries, max_steps ), phlist )
end
