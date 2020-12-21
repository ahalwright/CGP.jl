

function neutral_evolution( c::Chromosome, g::Goal, max_steps::Integer; print_steps::Bool=false )
  funcs = default_funcs( c.params.numinputs )
  step = 0
  ov = output_values( c) 
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  while step < max_steps && ov != g
    step += 1
    (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    new_ov = output_values( new_c )
    new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    #println("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur dis: ",current_distance,"  new_dis: ",new_distance )
    if new_ov == ov 
      c = new_c
      if print_steps
        println("step: ",step," is neutral.")
      end
    elseif new_distance < current_distance
      if print_steps
        println("step: ",step,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance)
      end
      c = new_c
      ov = new_ov
      current_distance = new_distance
    else
      if print_steps
        print("step: ",step,"  new_output: ",new_ov,"  new circuit: ")
        print_circuit( new_c )
      end 
    end
  end
  if step == max_steps
    println("neutral evolution failed with ",step," steps for goal: ",g)
    return nothing
  else
    println("neutral evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return c
  end
end

# Computes phenotype evolvability, genotype evolvability, robustness, complexity, steps
# num_circuits is the number of circuits used to compute properties.
function geno_circuits( g::Goal, p::Parameters, num_circuits::Integer, max_steps::Integer, max_attempts::Integer )
  funcs = default_funcs( p.numinputs )
  c = random_chromosome( p, funcs )
  circuit_list = Chromosome[]
  n_circuits = 0
  attempt = 0
  nc = nothing
  while attempt < max_attempts && n_circuits < num_circuits
    attempt += 1
    nc = neutral_evolution( c, g, max_steps )
    if nc != nothing
      n_circuits += 1
      push!( circuit_list, nc )
    end
  end
  if n_circuits < num_circuits
    println("geno_properties failed to find num_circuits circuits mapping to goal: ",g," in ", attempt," attempts.")
    return nothing
  end
  return circuit_list
end
 
function geno_properties( cl::Vector{Chromosome} )
  funcs = default_funcs( cl[1].params.numinputs )
  g = output_values(cl[1])
  sum_robust = 0.0
  sum_complexity = 0.0
  sum_evolvability = 0.0  # genotype evolability
  genotypes = Goal[]
  genotype_set = Set(Goal[])
  n_circuits = length(cl)
  for c in cl
    @assert output_values(c) == g
    sum_robust += mutational_robustness( c, funcs )
    sum_complexity += complexity5(c) 
    genotypes = Set(mutate_all(c,funcs,output_outputs=true,output_chromosomes=false))
    sum_evolvability += length(genotypes)-1
    #println("len(genotypes): ",length(genotypes),"  genotypes: ",genotypes)
    genotype_set = union( genotype_set, genotypes )
  end
  (g, sum_robust/n_circuits, sum_complexity/n_circuits, sum_evolvability/n_circuits, length(genotype_set) )
end

function geno_list_properties( gl::GoalList, p::Parameters, num_circuits::Integer, max_steps::Integer, max_attempts::Integer;
    csvfile::String="" )
  df = DataFrame()
  df.goal = Goal[]
  df.robustness = Float64[]
  df.complexity = Float64[]
  df.g_evolvability = Float64[]
  df.p_evolvability = Float64[]
  gp_list = Vector{Tuple{Goal,Float64,Float64,Float64,Float64}}[]
  gp_list = pmap( g->geno_properties( geno_circuits( g, p, num_circuits, max_steps, max_attempts ) ), gl )
  for gp in gp_list
    push!(df, gp) 
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# num_circuits: ",num_circuits)
      println(f,"# max_steps: ",max_steps)
      println(f,"# max_attempts: ",max_attempts)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
