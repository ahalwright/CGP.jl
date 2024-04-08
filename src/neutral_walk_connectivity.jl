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
# Note that sometimes the function fails when the neutral walk cannot be extended.
# There is another neutral_walk() function with the same signature in neutral_walk.jl
function neutral_walk( g::Goal, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  @assert p.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  circuit_ints = Int64[]
  n_repeats = 10  # maximum number of tries to find c that maps to g
  complexity_sum = 0.0
  complexity_count = 0
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  #println("complexity(c): ",complexity5(c))
  complexity_sum += complexity5(c)
  complexity_count += 1
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
      all_chromes = map( x->x[2], mutate_all( deepcopy(c), funcs, output_circuits=true )) 
      filter!( x->output_values(x)==g, all_chromes )  # only chromosomes that map to g
      #println(" attempts == maxtries  len(all_chromes): ",length(all_chromes))
      if length(all_chromes) == 0  
        println("unable to extend random_neutral_walk in function neutral_walk() at step: ",i)
        break
      else
        @assert output_values(all_chromes[1]) == g
        c = rand(all_chromes)
      end
    else
      c = new_c   # neutral walk extended
    end
    complexity_sum += complexity5(c)
    complexity_count += 1
    @assert output_values(c) == g
    # mutate_all() returns a list of pairs where the second element of the pair is the chromosome
    all_chromes = map( x->x[2], mutate_all( c, funcs, output_circuits=true )) 
    filter!( x->output_values(c)==g, all_chromes )
    for ch in all_chromes 
      push!( circuit_ints, circuit_int(ch) )
    end
    circuit_ints = unique(circuit_ints)  
  end # for loop
  (Set(circuit_ints), complexity_sum/complexity_count )
end
    
# Does multiple random walks and combines the returned circuit code lists if they have circuits in common.
function run_neutral_walk( g::Goal, p::Parameters, n_walks::Int64, steps::Int64, maxsteps::Int64, maxtries::Int64, int_list_file::IOStream )
  walk_list = Int64[]
  circuit_int_list = Set{Int64}[]
  walk_results = pmap(x->neutral_walk( g, p, steps, maxsteps, maxtries), collect(1:n_walks))  
  println("after pmap()")
  #walk_results = map(x->neutral_walk( g, p, steps, maxsteps, maxtries), collect(1:n_walks))  
  #println("walk_results[1]: ",walk_results[1])
  complexity_avg = 0.0
  for w = 1:n_walks
    #cclist = Set(neutral_walk( g, p, steps, maxsteps, maxtries))
    cclist = walk_results[w][1]
    complexity_avg += walk_results[w][2]
    #println("length(cclist): ",length(cclist))
    to_combine = Int64[]   # indices of circuit_int list to combine
    for i = 1:length(circuit_int_list)
      if length(intersect(cclist,circuit_int_list[i])) > 0
        push!(to_combine,i)
      end  
    end
    if length(to_combine) > 0
      println("w: ",w,"  to_combine: ",to_combine)
      combined_list = cclist
      for i = length(to_combine):-1:1
        #println("combining circuit_int_list[",to_combine[i],"]")
        combined_list = union(combined_list,circuit_int_list[to_combine[i]])
        if i > 1
          deleteat!( circuit_int_list, to_combine[i] )
          #println("removing index: ",to_combine[i]," from circuit_int_list")
          #println("walk_list before: ",walk_list)
          #walk_list = setdiff(walk_list,to_combine[i]) 
          deleteat!( walk_list, to_combine[i] )
          #println("removing index: ",to_combine[i]," from walk list")
          #println("walk_list after: ",walk_list)
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
    println("w: ",w,"  len: ",length(circuit_int_list),"  lengths circuit_int_list: ", [length(circuit_int_list[k]) for k = 1:length(circuit_int_list)])
    if w == n_walks
      println(int_list_file,"goal: ",g,"  len: ",length(circuit_int_list),"  lengths circuit_int_list: ", [length(circuit_int_list[k]) for k = 1:length(circuit_int_list)])
    end
    #print("w: ",w,"  length(circuit_int_list) after add: ",length(circuit_int_list))
    #println("  lengths circuit_int_list: ",[length(ccl) for ccl in circuit_int_list])
    for i = 1:length(circuit_int_list)
      for j = 1:(i-1)
        if !isempty(intersect(circuit_int_list[i],circuit_int_list[j]))
          println("XXXXX i: ",i,"  j: ",j,"  len intersect: ",length(intersect(circuit_int_list[i],circuit_int_list[j])))
        end
      end
    end
  end  # for w = 1:n_walks
  #walk_list
  (g,length(circuit_int_list),[length(ccl) for ccl in circuit_int_list],complexity_avg/n_walks)
end

# Do multiple runs of run_neutral_walk() defined above for each goal in goallist gl.
# To do repeated runs on a goal, include it multiple times in the gl.
function run_neutral_walk( gl::GoalList, p::Parameters, n_walks::Int64, steps::Int64, maxsteps::Int64, maxtries::Int64;
      csvfile::String="" )
  df = DataFrame()
  df.goal = Vector{MyInt}[]
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numints = Int64[]
  df.numlevsback = Int64[]
  df.n_walks = Int64[]
  df.steps = Int64[]
  df.maxsteps = Int64[]
  df.maxtries = Int64[]       
  df.n_combined = Float64[] 
  df.complexity = Float64[]
  if length(csvfile) > 0
    println("file: ","$(csvfile[1:(end-4)])_ints.txt")
    int_list_file=open("$(csvfile[1:(end-4)])_ints.txt","w")
  end
  result = map(g->run_neutral_walk( g, p, n_walks, steps, maxsteps, maxtries, int_list_file ), gl )
  #println("after map()")
  for r in result
    push!(df,(r[1], p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, n_walks, steps, maxsteps, maxtries, r[2], r[4] ))
  end
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  println("# date and time: ",Dates.now())
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )
  println("# funcs: ", Main.CGP.default_funcs(p.numinputs))
  print_parameters(p,comment=true)
  println("# n_walks: ",n_walks)
  println("# steps: ",steps)
  println("# maxsteps: ",maxsteps)
  println("# maxtries: ",maxtries)
  println("# ngoals: ",length(gl))
  if length(csvfile) > 0
    close(int_list_file)
    open( csvfile, "w" ) do f     
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# steps: ",steps)
      println(f,"# maxsteps: ",maxsteps)
      println(f,"# maxtries: ",maxtries)
      println(f,"# ngoals: ",length(gl))
      CSV.write(f, df, append=true, writeheader=true ) 
    end
  end
  return df
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
      all_chromes = map( x->x[2], mutate_all( deepcopy(c), funcs, output_circuits=true ))
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
    all_chromes = map( x->x[2], mutate_all( c, funcs, output_circuits=true ))
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

function neutral_walk_connectivity( p::Parameters,  g::Goal, max_steps::Int64, max_evolve_steps::Int64; c2_walk_length::Int64=200 )
  funcs = default_funcs( p.numinputs )
  res1 = mut_evolve( random_chromosome( p, funcs ), [g], funcs, max_evolve_steps )
  if res1[2] < max_evolve_steps
    c1 = res1[1]
  else
    error("failed to find circuit that maps to goal: ",g)
  end
  #print_build_chromosome( c1 )
  res2 = mut_evolve( random_chromosome( p, funcs ), [g], funcs, max_evolve_steps )
  if res2[2] < max_evolve_steps
    c2 = res2[1]
  else
    error("failed to find circuit that maps to goal: ",g)
  end
  #print_build_chromosome( c2 )
  neutral_walk_connectivity( c1, c2, c2_walk_length, max_steps )
end

function neutral_walk_connectivity( c1::Chromosome, c2::Chromosome, c2_walk_length::Int64, max_steps::Int64 )
  max_attempts = 10
  @assert c1.params.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  c2_walk_list = random_neutral_walk( c2, c2_walk_length )
  g = output_values(c1)
  out2 = output_values(c2)
  @assert g==out2
  c_code = circuit_code(c1)
  c_dist = circuit_distance_to_list( c1, c2_walk_list )
  println("original c_dist: ",c_dist)
  if c_dist == 0
    return 0
  end
  c = deepcopy(c1)
  @assert output_values(c) == g    # Assumes p.numoutputs==1
  for step = 1:max_steps
    c_code = circuit_code(c)
    (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    outputs = output_values(new_c)
    attempts = 1
    # try to find a c to extend walk
    while attempts < max_attempts && (outputs != g || circuit_distance_to_list(new_c,c2_walk_list) > c_dist || circuit_code(new_c)==c_code ) 
      (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      outputs = output_values(new_c)
      #println("step: ",step,"  attempts: ",attempts,"  outputs: ",outputs,"  new_c_dist: ",circuit_distance(new_c,c2))
      attempts += 1
    end
    if attempts == max_attempts
      # Alternatively, try all mutations of c in our attempt to extend walk
      c_code = circuit_code(c)
      all_chromes = map( x->x[2], mutate_all( deepcopy(c), funcs, output_circuits=true ))
      # filter to only chromosomes that map to g and do not increase circuit distance and change the circuit code
      filter!( (x->output_values(x)==g && circuit_distance_to_list(x,c2_walk_list) <= c_dist && circuit_code(x)!=c_code), all_chromes )  
      println(" attempts == maxtries  len(all_chromes): ",length(all_chromes))
      if length(all_chromes) == 0
        println("failed to extend random_neutral_walk in function neutral_walk() at step: ",step)
        return -1
      else
        @assert output_values(all_chromes[1]) == g
        c = rand(all_chromes)
        c_dist = circuit_distance_to_list(c,c2_walk_list) 
        println("mx step: ",step,"  neutral walk extended new_cdist: ",c_dist,"  new_circuit_code: ",circuit_code(c))
        @assert circuit_code(c) != c_code
      end
    else
      c = new_c   # neutral walk extended
      c_dist = circuit_distance_to_list(c,c2_walk_list) 
      println("at step: ",step,"  neutral walk extended new_cdist: ",c_dist,"  new_circuit_code: ",circuit_code(c))
      @assert circuit_code(c) != c_code
    end
    @assert output_values(c) == g
    new_c_dist = circuit_distance_to_list( c, c2_walk_list )
    @assert new_c_dist <= c_dist
    c_dist = new_c_dist < c_dist ? new_c_dist : c_dist
    if c_dist == 0.0
      println("found path c1 to c2 after ",step," steps")
      return step
    end
  end
  println("failed to find path from c1 to c2")
  return -1
end

function random_neutral_walk( c::Chromosome, length::Int64, max_attempts::Int64=100 )
  funcs = default_funcs( c.params.numinputs )
  chrome_list = Chromosome[c]
  g = output_values(c)
  c_code = circuit_code(c)
  new_c = Chromosome(p,[InputNode(1)],[],[OutputNode(1)],0.0,0.0)  # dummy chromosome to establish scope
  outputs = output_values(new_c)
  len = 0
  while len < length 
    attempts = 0
    while attempts < max_attempts && (outputs != g || circuit_code(new_c) == c_code )
      (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      outputs = output_values(new_c)
      #println("step: ",step,"  attempts: ",attempts,"  outputs: ",outputs )
      attempts += 1
    end
    if attempts < max_attempts
      c = new_c
      c_code = circuit_code(c)
      len += 1
      push!( chrome_list, new_c )
    else
      break
    end
  end
  return map( c->circuit_code(c), chrome_list )
end

function circuit_distance_to_list( c::Chromosome, c_code_list::Vector{Vector{Int64}} )
  circuit_distance_to_list( circuit_code(c), c_code_list )
end

function circuit_distance_to_list( c_code::Vector{Int64}, c_code_list::Vector{Vector{Int64}} )
  distance = length(c_code)
  for cc in c_code_list 
    @assert length(c_code) == length(cc)
    diff_count = 0
    for i = 1:length(c_code)
      diff_count += c_code[i] == cc[i] ? 0 : 1
    end
    distance = diff_count < distance ? diff_count : distance
  end
  distance/length(c_code)
end

# Added started in 12/17/21.  See notes/12_17_21.txt
# Use enumerate_circuits() in Chromosome.jl to get the list of Chromosomes (circuits)
function find_neutral_components( ch_list::Vector{Chromosome} )
  p = ch_list[1].params
  funcs = default_funcs(p.numinputs)
  S = Dict{Int128,Set{Int128}}()
  for g in ch_list
    S[chromosome_to_int(g)] = Set([chromosome_to_int(g)])
  end
  for g in ch_list
    ig = chromosome_to_int(g)
    ihlist = map(h->chromosome_to_int(h),mutate_all( g, funcs ) )
    for ih in ihlist
      if ih != ig
        S[ig] = union(S[ig],S[ih])
        s[ih] = S[ig]
      end
    end
  end
  S
end

