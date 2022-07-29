using DataFrames
using CSV
using Dates

# Compute the number of unique phenotypes in neighborhoods of increasing mutational size of a given genotype.
# Hypothesis:  Number unique phenotypes is postively related to complexity and 
#   negatively related to robustness.

# Parallel map run_mutational_evolvability() over the goals of gl.  Output a dataframe csv file.
# Note that goals in gl may be replaced by random goals if no circuit is evolved to map to a goal.
function run_mutational_evolvability( gl::GoalList, p::Parameters, num_mutations::Integer, max_steps::Integer,  max_attempts::Integer;
    csvfile::String="" )
  gl_sav = deepcopy(gl)
  println("gl: ",gl)
  #complexity_robustness_geno_counts = map( g->mutational_evolvability( g, p, num_mutations, max_steps, max_attempts ), gl )
  complexity_robustness_geno_counts = pmap( g->mutational_evolvability( g, p, num_mutations, max_steps, max_attempts ), gl )
  println("length pairs: ",length(complexity_robustness_geno_counts))
  df = DataFrame()
  df.goal = map(x->x[1],complexity_robustness_geno_counts)
  df.complexity = map(x->x[2],complexity_robustness_geno_counts)
  df.robusness = map(x->x[3],complexity_robustness_geno_counts)
  for i = 1:num_mutations
    insertcols!(df,size(df)[2]+1,Symbol("geno_count$i")=>map(x->x[i+3],complexity_robustness_geno_counts))
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# num_mutations: ",num_mutations)
      println(f,"# max_steps: ",max_steps)
      println(f,"# max_attempts: ",max_attempts)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
    
# Run the below defined mutational_evolvability() for a single goal.  
# Use neutral_evolution() to find the genotype that maps to the given goal (phenotype).
# If no circuit can be evolved that maps to g, then g is replace be a random goal (phenotype).
function mutational_evolvability( g::Goal, p::Parameters, num_mutations::Integer, max_steps::Integer,  max_attempts::Integer;
      max_goal_attempts::Int64=10)
  println("mut evolve g: ",g)
  nc = nothing
  goal_attempt = 0
  while goal_attempt < max_goal_attempts
    attempt = 0
    steps = max_steps
    while attempt < max_attempts && steps == max_steps
      c = random_chromosome(p,default_funcs(p.numinputs))
      attempt += 1
      (nc,steps) = neutral_evolution( c, funcs, g, max_steps )
      println("after neutral evolution: steps: ",steps)
    end
    if steps == max_steps
      g = randgoal(p.numinputs,p.numoutputs)
      println("new goal generated: ",g)
      goal_attempt += 1
    else 
      break
    end
  end
  complexity = complexity5(nc)
  #println("output values(nc): ",output_values(nc),"  complexity: ",complexity)
  print_circuit(nc)
  (robustness,goal_sets) = mutational_evolvability(nc, num_mutations )
  Tuple(vcat([g,complexity,robustness],[length(goal_sets[end-j]) for j = (num_mutations-1):-1:0]))
end
    
function mutational_evolvability( c::Chromosome, num_mutations::Int64 )
  funcs = default_funcs( c.params.numinputs )
  robustness = -1.0   # establish scope
  g = output_values(c)
  goal_sets = Set{Goal}[]
  circuit_set = Set{Chromosome}[]
  goal_set = Set( Goal[] )
  for i = 1:num_mutations
    if i == 1
      goal_circuit_pairs = mutate_all( c, funcs, output_circuits=true ) 
      circuit_set = Set( map( x->x[2], goal_circuit_pairs ))
      goal_list =  map( x->x[1], goal_circuit_pairs )
      robustness = length(filter(x->x==g,goal_list))/length(goal_list)
      goal_set = Set( goal_list )
    else
      circuit_set_previous = deepcopy( circuit_set )
      circuit_set = Set{Chromosome}[]
      goal_set = Set{Goal}[]
      j = 0
      for c in circuit_set_previous
        j += 1
        goal_circuit_pairs = mutate_all( c, funcs, output_circuits=true ) 
        circuit_set = union( Set(map( x->x[2], goal_circuit_pairs )), circuit_set )
        goal_set = union( Set( map( x->x[1], goal_circuit_pairs ) ), goal_set )
        #println("j: ",j,"  lengths: cl: ",length(circuit_list),"  gs: ",length(goal_set))
      end
    end    
    println("i: ",i,"  length(circuit_set): ",length(circuit_set))
    println("i: ",i,"  length(goal_set): ",length(goal_set))
    if i == 1
      push!(goal_sets,goal_set)
    else
      push!(goal_sets,union(goal_set,goal_sets[i-1]))
      println("i: ",i,"  goal: ",g,"  length(goal_sets[i]): ",length(goal_sets[i]))
    end
  end
  (robustness,goal_sets)
end
