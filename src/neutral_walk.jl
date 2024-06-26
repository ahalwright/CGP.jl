# Purpose: identify components of phenotype neutral network.

# steps is the number of random walk steps.
# maxsteps is the number of steps in mut_evolve().
# maxtries is the number of attempts to find the next chromosome in the random_neutral_walk.
# Assumes a single output goal
# Returns circuit_code_list of neutral circuits discovered on the neutral walk.
# Includes both walk circuits and neutral neighbors of walk circuits.
# There is another neutral_walk() function with the same signature in neutral_walk_connectivity.jl
function neutral_walk( g::Goal, p::Parameters, steps::Int64, maxsteps::Int64, maxtries::Int64 )
  n_repeats = 10  # maximum number of tries to find c that maps to g
  res = mut_evolve_repeat(n_repeats, p, [g], funcs, maxsteps )
  c = res[1]
  neutral_walk( c, p, steps, maxsteps, maxtries )
end

# neutral walk starting with a given chromosome instead of a given circuit.
# steps is the length of the neutral walk unless a walk isn't found
# maxtries is the number of mutations tried to find a neutral mutation
# The value of maxtries is not crucial because if a neutral mutation is not found in maxtries attempts,
#    all mutations are tried.
# maxsteps is never used so should be removed
function neutral_walk( c::Chromosome, steps::Int64, maxtries::Int64 )
  p = c.params
  g = output_values(c)
  @assert p.numoutputs == 1
  funcs = default_funcs( p.numinputs )
  walk_ints = Int64[]
  for i = 1:steps 
    (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    new_cint = chromosome_to_int( new_c )
    outputs = output_values(new_c)
    attempts = 1
    while attempts < maxtries && ((outputs != g) || !(new_cint in walk_ints))   # try to find a c to extend walk
      (new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      outputs = output_values(new_c)
      #println("attempts: ",attempts,"  outputs: ",outputs,"  goal: ",g)
      attempts += 1
    end
    if attempts == maxtries
      # Alternatively, try all mutations of c in our attempt to extend walk
      #println("alternatively")
      all_chromes = mutate_all( deepcopy(c), funcs, output_circuits=true, output_outputs=false )
      #println("all_chromes[1]: ",all_chromes[1])
      filter!( x->output_values(x)==g, all_chromes )  # only chromosomes that map to g
      #println(" attempts == maxtries  len(all_chromes): ",length(all_chromes))
      if length(all_chromes) == 0  
        error("unable to extend random_neutral_walk in function neutral_walk() at step: ",i)
      else
        new_c = rand(all_chromes)
        @assert output_values(new_c) == g
      end
    end
    @assert output_values(new_c) == g
    push!( walk_ints, chromosome_to_int(new_d) )
    #println("i: ",i,"  length(walk_ints): ",length(walk_ints),"  length(unique(walk_ints)): ",length(unique(walk_ints)))
  end # for loop
  unique(walk_ints)
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

# Extract the list of discovered component sizes from an *_ints.txt file created by neutral_walks() function in this file
# Return a dataframe with two fields:  goal, and comp_list.
function extract_component_lists( intsfile::String )
  csdf = DataFrame()
  csdf.goal = String[]
  csdf.comp_list = String[]
  open( intsfile, "r" ) do f
    for line in readlines(f)
      m = match(r"(\[0x.*\]).*(\[.*\])",line) 
      push!(csdf,m.captures)
    end
  end
  csdf
end

# Add the component sizes lists returned by function extract_component_sizes() to the dataframe of csvfile.
#  The resulting dataframe is written to a csvfile with "cs" appended to the file name (before the ".csv")
function add_component_sizes_to_df( csvfile::String )
  df = read_dataframe(csvfile)
  ints_file = string( csvfile[1:(end-4)], "_ints.txt" )
  csdf = extract_component_lists( ints_file )
  size_lists = [ eval(Meta.parse(csdf.comp_list[i])) for i = 1:size(csdf)[1] ]
  size_lists = map(x->sort(x,rev=true),size_lists)
  df.size_sum = map( sum, size_lists )
  df.size_ratio = [ (df.size_sum[i] - size_lists[i][1])/df.size_sum[i] for i = 1:size(df)[1] ]
  df.comp_list= csdf.comp_list
  dfcs_file = string( csvfile[1:(end-4)], "ccs.csv" )
  println("dfcs_file: ",dfcs_file )
  write_dataframe_with_comments(df,csvfile,dfcs_file)
end
