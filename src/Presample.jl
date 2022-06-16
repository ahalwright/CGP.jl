# Goal:  Evolve genotypes to map to a given target phenotype where this is difficult, such as 5 inputs.
# Method:  Build a random pre-sample of genotype-phenotype pairs by sampling genotypes and computing their phenotypes.
# Choose a subset of the pre-sample whose phenotypes are close to the target phenotype, and start epochal evolutions from this subset.

# Do both psample_evolve() and run_neutral_evol() and save result to csvfile
function psample_neutral_evolve( p::Parameters, funcs::Vector{Func}, target_pheno::Goal, psample_size::Int64, num_evolves::Int64, max_steps::Int64;
    use_lincircuit::Bool=false, csvfile::String="" )
  (pse,ptime) = @timed psample_evolve( p, funcs, rph, psample_size, num_evolves, max_steps )
  (rlist,ntime) = @timed run_neutral_evol( p, funcs, rph, num_evolves, max_steps )
  if length(csvfile) > 0
    df = DataFrame()
    df.psample_evolve = map(x->x[2],pse)
    df.neutral_evolve = map(x->x[2],rlist) 
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# target phenotype: ",rph)
      println(f,"# num_evolves: ",num_evolves)
      println(f,"# psample_size: ",psample_size)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
      return df
    end
  end
end

function psample_evolve( p::Parameters, funcs::Vector{Func}, target_pheno::Goal, psample_size::Int64, num_evolves::Int64, max_steps::Int64;
    use_lincircuit::Bool=false )
  psample = pre_sample( p, funcs, psample_size, use_lincircuit=use_lincircuit )
  best_sub_sample = best_subsample( p, psample, target_pheno, num_evolves )
  solution_list = pmap( bss->neutral_evolution( bss[1][1], target_pheno, max_steps, funcs=funcs ), best_sub_sample )  # bss[1][1] is the circuit
  #solution_list = map( bss->neutral_evolution( bss[1][1], target_pheno, max_steps, funcs=funcs ), best_sub_sample )  # bss[1][1] is the circuit
end

function pre_sample( p::Parameters, funcs::Vector{Func}, sample_size::Int64; use_lincircuit::Bool=false )
  pre_sample_list = Tuple{Circuit,Goal}[]
  for i = 1:sample_size
    if use_lincircuit
      c = rand_lcircuit( p, funcs )
    else
      c = random_chromosome( p, funcs )
    end
    ph = output_values(c)
    push!(pre_sample_list,(c,ph))
  end
  pre_sample_list
end
  
function best_subsample( p::Parameters, pre_sample_list::Vector{Tuple{Circuit,Goal}}, target::Goal, subsample_size::Int64 )
  pre_sample_with_hamming_distances = map( s->(s,hamming_distance( target, s[2], p.numinputs )), pre_sample_list )
  sort!( pre_sample_with_hamming_distances, by=x->x[2], rev=true )
  pre_sample_with_hamming_distances[1:subsample_size]
end

# Note that Ones is a global variable set by calling default_funcs(p)
function rand_pheno( p::Parameters )
  [ rand(MyInt(0):Ones) for i = 1:p.numoutputs ]
end

function run_neutral_evol( p::Parameters, funcs::Vector{Func}, target_pheno::Goal, num_evolves::Int64, max_steps::Int64;
    use_lincircuit::Bool=false )
  result_list = pmap( i->neutral_evolution( (use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs ) ), target_pheno, max_steps, funcs=funcs ), 1:num_evolves )
end
