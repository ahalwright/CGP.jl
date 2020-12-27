using DataFrames
using CSV
using Dates

function mutational_evolvability( c::Chromosome, num_mutations::Int64 )
  funcs = default_funcs( c.params.numinputs )
  goal_sets = Set{Goal}[]
  circuit_set = Set{Chromosome}[]
  goal_set = Set( Goal[] )
  for i = 1:num_mutations
    if i == 1
      goal_circuit_pairs = mutate_all( c, funcs, output_chromosomes=true ) 
      circuit_set = Set( map( x->x[2], goal_circuit_pairs ))
      goal_set = Set( map( x->x[1], goal_circuit_pairs ) )
    else
      circuit_set_previous = deepcopy( circuit_set )
      circuit_set = Set{Chromosome}[]
      goal_set = Set{Goal}[]
      j = 0
      for c in circuit_set_previous
        j += 1
        goal_circuit_pairs = mutate_all( c, funcs, output_chromosomes=true ) 
        circuit_set = union( Set(map( x->x[2], goal_circuit_pairs )), circuit_set )
        goal_set = union( Set( map( x->x[1], goal_circuit_pairs ) ), goal_set )
        #println("j: ",j,"  lengths: cl: ",length(circuit_list),"  gs: ",length(goal_set))
      end
    end    
    println("length(circuit_set): ",length(circuit_set))
    println("length(goal_set): ",length(goal_set))
    if i == 1
      push!(goal_sets,goal_set)
    else
      push!(goal_sets,union(goal_set,goal_sets[i-1]))
      println("i: ",i,"  length(goal_sets[i]): ",length(goal_sets[i]))
    end
  end
  goal_sets 
end
    
function run_mutational_evolvability( g::Goal, p::Parameters, num_mutations::Integer, max_steps::Integer,  max_attempts::Integer)
  nc = nothing
  attempt = 0
  steps = 0
  c = random_chromosome(p,default_funcs(p.numinputs))
  while attempt < max_attempts && steps == max_steps
    attempt += 1
    (nc,steps) = neutral_evolution( c, g, max_steps, max_attempts )
  end
  gs=mutational_evolvability(c, num_mutations )
  length(gs[end])
end
    
function run_mutational_evolvability( gl::GoalList, p::Parameters, num_mutations::Integer, max_steps::Integer,  max_attempts::Integer;
    csvfile::String="" )
  geno_counts = pmap( g->run_mutational_evolvability( g, p, num_mutations, max_steps, max_attempts ), gl )
  df = DataFrame()
  df.goal = gl
  df.geno_counts = geno_counts
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
  df
end
