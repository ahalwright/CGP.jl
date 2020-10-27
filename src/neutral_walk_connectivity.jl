# Objective:  Find the size of the components of the neutral network of a goal.
# Methodology:  Do multiple neutral random walks starting at a random circuit that maps to the goal.
# neutral_walk() accumulate all neighbors of all neutral circuits encountered on the random neutral walk.
# Then run_neutral_walk() combines the circuit lists returned by multiple random walks 
#      if they have a circuit in common.
# Thus, run_neutral_walk() returns circuit lists that have not been shown to be in the same connected
#      component.

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
    sav_c = deepcopy(c)
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
