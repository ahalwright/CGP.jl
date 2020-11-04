# Objective:  Find the size of the components of the neutral network of a goal.
# Methodology:  Do multiple neutral random walks starting at a random circuit that maps to the goal.
# neutral_walk() accumulate all neighbors of all neutral circuits encountered on the random neutral walk.
# Then run_neutral_walk() combines the circuit lists returned by multiple random walks 
#      if they have a circuit in common.
# Thus, run_neutral_walk() returns circuit lists that have not been shown to be in the same connected
#      component.

using Statistics, Printf

# steps is the number of random walk steps.
# maxsteps is the number of steps in mut_evolve().
# maxtries is the number of attempts to find the next chromosome in the random_neutral_walk.
# Assumes a single output goal
# Returns circuit_code_list of neutral circuits discovered on the neutral walk.
# Includes both walk circuits and neutral neighbors of walk circuits.
function neutral_walk( g::Goal, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  @assert p.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  circuit_codes = Vector{Int64}[]
  n_repeats = 10  # maximum number of tries to find c that maps to g
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  println("complexity(c): ",complexity5(c))
  #@assert sort(output_values(c)) == sort(g)
  @assert output_values(c) == g    # Assumes p.numoutputs==1
  push!( circuit_codes, chromosome_code(c) )
  for i = 1:steps 
    (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    outputs = output_values(new_c)
    attempts = 1
    while attempts < maxtries && outputs != g  # try to find a c to extend walk
      (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      outputs = output_values(new_c)
      #println("attempts: ",attempts,"  outputs: ",outputs,"  goal: ",g)
      attempts += 1
    end
    if attempts == maxtries
      # Alternatively, try all mutations of c in our attempt to extend walk
      all_chromes = map( x->x[2], mutate_all( deepcopy(c), funcs, output_chromosomes=true )) 
      filter!( x->output_values(x)==g, all_chromes )  # only chromosomes that map to g
      #println(" attempts == maxtries  len(all_chromes): ",length(all_chromes))
      if length(all_chromes) == 0  
        println("unable to extend random_neutral_walk in function neutral_walk() at step: ",i)
        continue
      else
        @assert output_values(all_chromes[1]) == g
        c = rand(all_chromes)
      end
    else
      c = new_c   # neutral walk extended
    end
    @assert output_values(c) == g
    # mutate_all() returns a list of pairs where the second element of the pair is the chromosome
    all_chromes = map( x->x[2], mutate_all( c, funcs, output_chromosomes=true )) 
    filter!( x->output_values(c)==g, all_chromes )
    for ch in all_chromes 
      push!( circuit_codes, chromosome_code(ch) )
    end
    circuit_codes = unique(circuit_codes)  
  end # for loop
  circuit_codes
end
    
# Does multiple random walks and combines the returned circuit code lists if they have circuits in common.
function run_neutral_walk( g::Goal, p::Parameters, n_walks::Int64, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  walk_list = Int64[]
  circuit_code_list = Vector{Vector{Int64}}[]
  for w = 1:n_walks
    cclist = neutral_walk( g, p, steps, maxsteps, maxtries)
    println("length(cclist): ",length(cclist))
    to_combine = Int64[]   # indices of circuit_code list to combine
    for i = 1:length(circuit_code_list)
      if length(intersect(cclist,circuit_code_list[i])) > 0
        #new_list = unique(union(cclist,cclist))
        push!(to_combine,i)
      end  
    end
    println("to_combine: ",to_combine)
    if length(to_combine) > 0
      combined_list = cclist
      for i = length(to_combine):-1:1
        println("combining circuit_code_list[",to_combine[i],"]")
        combined_list = unique(union(combined_list,circuit_code_list[to_combine[i]]))
        if i > 1
          deleteat!( circuit_code_list, to_combine[i] )
          println("removing index: ",to_combine[i]," from circuit_code_list")
          println("walk_list before: ",walk_list)
          #walk_list = setdiff(walk_list,to_combine[i]) 
          deleteat!( walk_list, to_combine[i] )
          println("removing index: ",to_combine[i]," from walk list")
          println("walk_list after: ",walk_list)
        end
      end
      circuit_code_list[to_combine[1]] = combined_list
      #=
      for j = 2:length(to_combine)
        deleteat!( circuit_code_list, j )
      end
      println("walk_list before: ",walk_list)
      walk_list = setdiff(walk_list,to_combine[2:end])
      println("walk_list after: ",walk_list)
      =#
    else
      push!(walk_list,w)
      push!(circuit_code_list,cclist)
    end
    println("w: ",w,"  walk_list after add: ",walk_list)
    #println("  lengths circuit_code_list: ",[length(circuit_code_list[walk_list[j]]) for j = 1:length(walk_list)])
    print("w: ",w,"  length(circuit_code_list) after add: ",length(circuit_code_list))
    println("  lengths circuit_code_list: ",[length(ccl) for ccl in circuit_code_list])
    for i = 1:length(circuit_code_list)
      for j = 1:(i-1)
        if !isempty(intersect(circuit_code_list[i],circuit_code_list[j]))
          println("XXXXX i: ",i,"  j: ",j,"  len intersect: ",length(intersect(circuit_code_list[i],circuit_code_list[j])))
        end
      end
    end
  end
  walk_list
end

# Distribution of complexities on a neutral walk.   
function neutral_walk_complexity( g::Goal, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  @assert p.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  walk_complexities = Float64[]
  walk_count = 0
  neighbor_complexities = Float64[]
  neighbor_count = 9
  n_repeats = 10  # maximum number of tries to find c that maps to g
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  complexity = complexity5(c)
  #println("w complexity(c): ",complexity)
  push!(walk_complexities,complexity)
  #@assert sort(output_values(c)) == sort(g)
  @assert output_values(c) == g    # Assumes p.numoutputs==1
  for i = 1:steps
    (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    outputs = output_values(new_c)
    attempts = 1
    while attempts < maxtries && outputs != g  # try to find a c to extend walk
      (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      outputs = output_values(new_c)
      #println("attempts: ",attempts,"  outputs: ",outputs,"  goal: ",g)
      attempts += 1
    end
    if attempts == maxtries
      # Alternatively, try all mutations of c in our attempt to extend walk
      all_chromes = map( x->x[2], mutate_all( deepcopy(c), funcs, output_chromosomes=true ))
      filter!( x->output_values(x)==g, all_chromes )  # only chromosomes that map to g
      #println(" attempts == maxtries  len(all_chromes): ",length(all_chromes))
      if length(all_chromes) == 0
        println("unable to extend random_neutral_walk in function neutral_walk() at step: ",i)
        continue
      else
        @assert output_values(all_chromes[1]) == g
        c = rand(all_chromes)
      end
    else
      c = new_c   # neutral walk extended
    end     
    @assert output_values(c) == g
    complexity = complexity5(c)
    #println("w complexity(c): ",complexity)
    push!(walk_complexities,complexity)
    walk_count += 1
    # mutate_all() returns a list of pairs where the second element of the pair is the chromosome
    all_chromes = map( x->x[2], mutate_all( c, funcs, output_chromosomes=true ))
    filter!( x->output_values(c)==g, all_chromes )
    for ch in all_chromes
      complexity = complexity5(ch)
      #println("n complexity(c): ",complexity)
      push!(neighbor_complexities,complexity)
      neighbor_count += 1
    end
  end
  #println("steps: ",steps,"  length(walk_c): ",length(walk_complexities),"  length(nbr_c): ",length(neighbor_complexities))
  (walk_complexities,neighbor_complexities)
end

# Distribution of complexities on a neutral walk.   
function run_neutral_walk_complexity( goallist::GoalList, p::Parameters, n_walks::Int64, steps::Int64, 
      maxsteps::Int64, maxtries::Int64; csvfile::String="" )
  if length(csvfile) > 0
    f = open(csvfile,"w")
    print_parameters(f,p)
    println(f,"n_walks: ",n_walks)
    println(f,"steps: ",steps)
    println(f,"maxtries: ",maxtries)
  end
  for g in goallist
    w_cmplx = fill(Float64[],n_walks)
    n_cmplx = fill(Float64[],n_walks)
    walk_complexities = Float64[]
    neighbor_complexities = Float64[]
    for w = 1:n_walks
      (w_cmplx[w],n_cmplx[w]) = neutral_walk_complexity( g, p, steps, maxsteps, maxtries )
      #println("w: ",w,"  length(w_cmplx): ",length(w_cmplx[w]),"  length(n_cmplx): ",length(n_cmplx[w]))
    end
    walk_complexities = vcat( w_cmplx... )
    neighbor_complexities = vcat( n_cmplx... )
    println("length(walk_complexities): ",length(walk_complexities),"  length(neighbor_complexities): ",length(neighbor_complexities))
    # Statistics
    w_mean = mean(walk_complexities) 
    w_max = maximum(walk_complexities)     
    w_min = minimum(walk_complexities)
    w_std = std(walk_complexities)
    w_q = quantile(walk_complexities,[0.9,0.95,0.99])
    n_mean = mean(neighbor_complexities) 
    n_max = maximum(neighbor_complexities)
    n_min = minimum(neighbor_complexities)     
    n_std = std(neighbor_complexities)
    n_q = quantile(neighbor_complexities,[0.9,0.95,0.99])
    @printf("goal: [0x%04x]",g[1])
    @printf("  w_count:%8d",length(walk_complexities))
    @printf("  w_mean:%6.3f",w_mean)
    @printf("  w_max:%6.3f",w_max)
    @printf("  w_min:%6.3f",w_min)
    @printf("  w_std:%6.3f",w_std)
    @printf("  w_q90:%6.3f",w_q[1])
    @printf("  w_q95:%6.3f",w_q[2])
    @printf("  w_q99:%6.3f\n",w_q[3])
    @printf("goal: [0x%04x]",g[1])
    @printf("  n_count:%8d",length(neighbor_complexities))
    @printf("  n_mean:%6.3f",n_mean)
    @printf("  n_max:%6.3f",n_max)
    @printf("  n_min:%6.3f",n_min)
    @printf("  n_std:%6.3f",n_std)
    @printf("  n_q90:%6.3f",n_q[1])
    @printf("  n_q95:%6.3f",n_q[2])
    @printf("  n_q99:%6.3f\n",n_q[3])
    if length(csvfile) > 0
      @printf(f,"goal: [0x%04x]",g[1])
      @printf(f,"  w_count:%8d",length(walk_complexities))
      @printf(f,"  w_mean:%6.3f",w_mean)
      @printf(f,"  w_max:%6.3f",w_max)
      @printf(f,"  w_min:%6.3f",w_min)
      @printf(f,"  w_std:%6.3f",w_std)
      @printf(f,"  w_q90:%6.3f",w_q[1])
      @printf(f,"  w_q95:%6.3f",w_q[2])
      @printf(f,"  w_q99:%6.3f\n",w_q[3])
      @printf(f,"goal: [0x%04x]",g[1])
      @printf(f,"  n_count:%8d",length(neighbor_complexities))
      @printf(f,"  n_mean:%6.3f",n_mean)
      @printf(f,"  n_max:%6.3f",n_max)
      @printf(f,"  n_min:%6.3f",n_min)
      @printf(f,"  n_std:%6.3f",n_std)
      @printf(f,"  n_q90:%6.3f",n_q[1])
      @printf(f,"  n_q95:%6.3f",n_q[2])
      @printf(f,"  n_q99:%6.3f\n",n_q[3])
    end
  end
  close(f)
  #(walk_complexities,neighbor_complexities)
end
    
