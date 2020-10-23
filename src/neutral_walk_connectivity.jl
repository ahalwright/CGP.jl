
# steps is the number of random walk steps.
# maxsteps is the number of steps in mut_evolve().
# maxtries is the number of attempts to find the next chromosome in the random_neutral_walk.
# Assumes a single output goal
function neutral_walk( g::Goal, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  @assert p.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  circuit_code_list = Vector{Int64}[]
  n_repeats = 10  # maximum number of tries to find c that maps to g
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  #@assert sort(output_values(c)) == sort(g)
  @assert output_values(c) == g    # Assumes p.numoutputs==1
  push!( circuit_code_list, chromosome_code(c) )
  for i = 1:steps 
    sav_c = deepcopy(c)
    (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    outputs = output_values(new_c)
    attempts = 1
    while attempts < maxtries && outputs != g 
      (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      outputs = output_values(new_c)
      #println("attempts: ",attempts,"  outputs: ",outputs,"  goal: ",g)
      attempts += 1
    end
    if attempts == maxtries
      all_chromes = map( x->x[2], mutate_all( deepcopy(c), funcs, output_chromosomes=true )) 
      filter!( x->output_values(x)==g, all_chromes )
      println(" attempts == maxtries  len(all_chromes): ",length(all_chromes))
      if length(all_chromes) == 0
        println("unable to extend random_neutral_walk in function neutral_walk() at step: ",i)
        continue
      else
        @assert output_values(all_chromes[1]) == g
        c = rand(all_chromes)
      end
    else
      c = new_c
    end
    @assert output_values(c) == g
    # mutate_all() returns a list of pairs where the second element of the pair is the chromosome
    all_chromes = map( x->x[2], mutate_all( c, funcs, output_chromosomes=true )) 
    filter!( x->output_values(c)==g, all_chromes )
    for ch in all_chromes 
      push!( circuit_code_list, chromosome_code(ch) )
    end
    circuit_code_list = unique(circuit_code_list)
  end
  circuit_code_list
end
    
function run_neutral_walk( g::Goal, p::Parameters, n_walks::Int64, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  walk_list = Int64[]
  circuit_code_list = Vector{Vector{Int64}}[]
  for w = 1:n_walks
    nwlist = neutral_walk( g, p, steps, maxsteps, maxtries)
    #push!(walk_list,w)
    #new_list = Vector{Vector{Int64}}[]
    to_combine = Int64[]
    for i = 1:length(circuit_code_list)
      if length(intersect(nwlist,circuit_code_list[i])) > 0
        #new_list = unique(union(nwlist,cclist))
        push!(to_combine,i)
      end  
    end
    println("to_combine: ",to_combine)
    if length(to_combine) > 0
      combined_list = nwlist
      for i = 1:length(to_combine)
        println("combining circuit_code_list[",i,"]")
        combined_list = unique(union(combined_list,circuit_code_list[i]))
      end
      circuit_code_list[to_combine[1]] = combined_list
      for j = 2:length(to_combine)
        deleteat!( circuit_code_list, j )
        println("removing index: ",j," from circuit_code_list")
      end
      println("walk_list before: ",walk_list)
      walk_list = setdiff(walk_list,to_combine[2:end])
      println("walk_list after: ",walk_list)
    else
      push!(walk_list,w)
      #println("walk_list after add: ",walk_list)
      push!(circuit_code_list,nwlist)
      #println("length(circuit_code_list) after add: ",length(circuit_code_list))
    end
    for i = 1:length(circuit_code_list)
      for j = 1:(i-1)
        if !isempty(intersect(circuit_code_list[i],circuit_code_list[j]))
          println("i: ",i,"  j: ",j,"len intersect: ",length(intersect(circuit_code_list[i],circuit_code_list[j])))
        end
      end
    end
  end
  walk_list
end
