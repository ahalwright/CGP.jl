export robustness, genotype_robustness, phenotype_robustness

function robustness( c::Circuit, funcs::Vector{Func} )
  #print("robustness: c:  ")
  #print_circuit(c,funcs)
  c_output = output_values(c,funcs)
  outputs = mutate_all( c, funcs, output_outputs=true )
  #println("outputs[1]: ",outputs[1])
  robust_outputs = filter( x->x==c_output, outputs )
  return length(robust_outputs)/length(outputs)
end   

# Example:  genotype_robustness( [(3,8,4),(4,10,5)], 100 )
function genotype_robustness( ni_ng_lb_triples::Vector{Tuple{Int64,Int64,Int64}}, nreps::Int64 )
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
    for i = 1:nreps
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
  df = DataFrame( :ninteriors=>ni_list, :ngates=>ng_list, :lb=>lb_list, :nreps=>fill(nreps,length(ni_ng_lb_triples)), 
      :mean_robust=>mean_rbst_list, :std_robust=>std_rbst_list, :q10_robust=>q10_rbst_list, :q90_robust=>q90_rbst_list,
      :rlists=>rlist_list )
end

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
  for (ni, ng, lb) in ni_ng_lb_triples
    p = Parameters( ni, 1, ng, lb )
    funcs = default_funcs( p )
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
  df = DataFrame( :ninteriors=>ni_list, :ngates=>ng_list, :lb=>lb_list, :numcircuits=>fill(numcircuits,length(ni_ng_lb_triples)), :ngoals=>fill(ngoals,length(ni_ng_lb_triples)), 
        :mean_robust=>mean_rbst_list, :std_robust=>std_rbst_list, :q10_robust=>q10_rbst_list, :q90_robust=>q90_rbst_list )
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
end

# Evolves numciruits circuits to map to each phenotype in phlist.
# Returns a dataframe containing robustness statistics and rlists.
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
