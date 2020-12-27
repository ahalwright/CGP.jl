using Statistics, Printf

# Calculates robustness and evolvability for all goal functions by sampling circuits as in Hu et al. (2020).
# The node-adjacency matrix for the network of phenotypes is computed as goal_edge_matrix.
# Diagonal entries are the number of self edges which is the robustness count.
# goal_edge_matrix[g,h] is the number of mutations discovered from goal g to goal h.
# The methodology is to do nprocesses*nwalks random walks each of length steps, and record all of the transitions made.
# If csvfile is specified, writes the dataframe to this file.
function run_random_walks_parallel( nprocesses::Int64, nwalks::Int64, p::Parameters, steps::Int64; csvfile::String="" )
  ngoals = 2^(2^p.numinputs)
  goal_edge_matrix = zeros(Int64,ngoals,ngoals)
  goal_edge_matrix_list = pmap( x->run_random_walks( nwalks, p, steps ), collect(1:nprocesses) )
  println("len: ",length(goal_edge_matrix_list),"  size: ",size(goal_edge_matrix_list[1]),"  type: ",typeof(goal_edge_matrix_list[1]))
  for gem in goal_edge_matrix_list
    goal_edge_matrix .+= gem
  end
  df = robust_evolvability( goal_edge_matrix )
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# nwalks: ",nwalks)
      println(f,"# steps: ",steps)
      println(f,"# nprocesses: ",nprocesses)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
end

function random_walk( c::Chromosome, steps::Int64 )
  funcs = default_funcs(c.params.numinputs)
  ngoals = 2^(2^c.params.numinputs)
  goal_edge_matrix = zeros(Int64,ngoals,ngoals)
  goal_pair_list = Tuple{MyInt,MyInt}[]
  goal = output_values(c)[1]
  #println("start goal: ",goal)
  for i = 1:steps 
    (new_c,active) = mutate_chromosome!( c, funcs )
    prev_goal = goal  
    goal = output_values(c)[1]
    #println("i: ",i,"  goal: ",goal)
    goal_edge_matrix[Int64(prev_goal)+1,Int64(goal)+1] += 1
    push!( goal_pair_list, (prev_goal,goal) )
  end
  goal_edge_matrix
end

function run_random_walks( nwalks::Int64, p::Parameters, steps::Int64 )
  println("run_random_walks with nwalks: ",nwalks,"  steps: ",steps)
  @assert p.numoutputs == 1
  funcs = default_funcs(p.numinputs)
  ngoals = 2^(2^p.numinputs)
  goal_edge_matrix = zeros(Int64,ngoals,ngoals)
  for i = 1:nwalks
    c = random_chromosome( p, funcs )
    goal_edge_matrix .+= random_walk( c, steps )
  end
  goal_edge_matrix 
end

# Calculates robustness and degree evolvability for each goal and saves these in a dataframe.
# Robustness 
function robust_evolvability( goal_edge_matrix::Array{Int64,2} )
  ngoals = size(goal_edge_matrix)[1]
  triangularize!(goal_edge_matrix)
  robustness = zeros(Float64,ngoals)
  evolvability = zeros(Float64,ngoals)
  for g = 1:ngoals
    robustness[g] = goal_edge_matrix[g,g]
    sum1 = g==ngoals ? 0 : sum( ((goal_edge_matrix[g,h] != 0) ? 1 : 0) for h = (g+1):ngoals )
    pairs1 = g==ngoals ? 0 : [ (g,h) for h = (g+1):ngoals ] 
    sum2 = g==1 ? 0 : sum( ((goal_edge_matrix[h,g] != 0) ? 1 : 0) for h = 1:(g-1))
    pairs2 = g==1 ? 0 : [ ((h,g),goal_edge_matrix[h,g]) for h = 1:(g-1)] 
    #println("g: ",g,"  pairs1: ",pairs1,"  sum1: ",sum1)
    #println("g: ",g,"  pairs2: ",pairs2,"  sum2: ",sum2)
    evolvability[g] = sum1 + sum2
  end
  df = DataFrame()
  df.goal = collect(MyInt(0):MyInt(ngoals-1))
  df.robustness = robustness
  df.evolvability = evolvability
  df
end

# Convert into a upper triangular matrix by adding the below-diagonal entries to the corresponding above-diagonal entry
function triangularize!( M::Array{Int64,2} )
  ngoals = size(M)[1]
  @assert ngoals == size(M)[2]
  for g = 1:(ngoals-1)
    for h = (g+1):ngoals
      M[g,h] += M[h,g]
      M[h,g] = 0
    end
  end   
end

# steps is the number of random walk steps.
# maxsteps is the number of steps in mut_evolve().
# maxtries is the number of attempts to find the next chromosome in the random_neutral_walk.
# Assumes a single output goal
# Returns circuit_code_list of neutral circuits discovered on the neutral walk.
# Includes both walk circuits and neutral neighbors of walk circuits.
function neutral_walk( g::Goal, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  n_repeats = 10  # maximum number of tries to find c that maps to g
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  neutral_walk( c, p, steps, maxsteps, maxtries )
end

# neutral walk starting with a given chromosome instead of a given circuit.
function neutral_walk( c::Chromosome, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  @assert p.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  circuit_ints = Int64[]
  n_repeats = 10  # maximum number of tries to find c that maps to g
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  println("complexity(c): ",complexity5(c))
  #@assert sort(output_values(c)) == sort(g)
  @assert output_values(c) == g    # Assumes p.numoutputs==1
  push!( circuit_ints, circuit_int(c) )
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
        error("unable to extend random_neutral_walk in function neutral_walk() at step: ",i)
        #continue
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
      push!( circuit_ints, circuit_int(ch) )
    end
    circuit_ints = unique(circuit_ints)  
  end # for loop
  circuit_ints
end
    
# Does multiple random walks and combines the returned circuit code lists if they have circuits in common.
function run_neutral_walk( g::Goal, p::Parameters, n_walks::Int64, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  walk_list = Int64[]
  circuit_int_list = Set{Int64}[]
  for w = 1:n_walks
    cclist = Set(neutral_walk( g, p, steps, maxsteps, maxtries))
    println("length(cclist): ",length(cclist))
    to_combine = Int64[]   # indices of circuit_int list to combine
    for i = 1:length(circuit_int_list)
      if length(intersect(cclist,circuit_int_list[i])) > 0
        push!(to_combine,i)
      end  
    end
    println("to_combine: ",to_combine)
    if length(to_combine) > 0
      combined_list = cclist
      for i = length(to_combine):-1:1
        println("combining circuit_int_list[",to_combine[i],"]")
        combined_list = union(combined_list,circuit_int_list[to_combine[i]])
        if i > 1
          deleteat!( circuit_int_list, to_combine[i] )
          println("removing index: ",to_combine[i]," from circuit_int_list")
          println("walk_list before: ",walk_list)
          #walk_list = setdiff(walk_list,to_combine[i]) 
          deleteat!( walk_list, to_combine[i] )
          println("removing index: ",to_combine[i]," from walk list")
          println("walk_list after: ",walk_list)
        end
      end
      circuit_int_list[to_combine[1]] = combined_list
    else
      push!(walk_list,w)
      push!(circuit_int_list,cclist)
      #println("walk list: ",walk_list)
      #println("   lengths circuit_int_list: ", [length(circuit_int_list[k]) for k = 1:length(circuit_int_list)])
    end
    #println("w: ",w,"  walk_list after add: ",walk_list)
    #println("   lengths circuit_int_list: ", [length(circuit_int_list[k]) for k = 1:length(circuit_int_list)])
    print("w: ",w,"  length(circuit_int_list) after add: ",length(circuit_int_list))
    println("  lengths circuit_int_list: ",[length(ccl) for ccl in circuit_int_list])
    for i = 1:length(circuit_int_list)
      for j = 1:(i-1)
        if !isempty(intersect(circuit_int_list[i],circuit_int_list[j]))
          println("XXXXX i: ",i,"  j: ",j,"  len intersect: ",length(intersect(circuit_int_list[i],circuit_int_list[j])))
        end
      end
    end
  end
  walk_list
end
