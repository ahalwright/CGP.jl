
function run_random_neutral_walk( p::Parameters, ngoals::Int64, steps::Int64, maxsteps::Int64 )
  df = DataFrame()
  df.goal = Goal[]
  df.evolvable_count = Int64[]
  df.complexity = Float64[]
  funcs = default_funcs( p.numinputs )
  for gc = 1:ngoals
    g = randgoal( p.numinputs, p.numoutputs )
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve_repeat( 10, p, [g], funcs, maxsteps)
    res = mut_evolve_repeat( 10, p, [g], funcs, maxsteps)
    if res == nothing
      error("run_random_neutral_walk: mut_evolve_repeat() failed.  Try increasing maxsteps.")
    end
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = res
    @assert output_values(c) == g
    complexity = complexity5(c)
    ev_count = random_neutral_walk( c, steps, maxsteps, evo_count_only=true )
    push!(df, (g, ev_count, complexity ))
  end
  df
end

# See diary10_13.txt for example runs
function run_random_neutral_walk( p::Parameters, g::Goal, steps::Int64, maxsteps::Int64, nreps::Int64=1; skip_steps::Int64=1 )
  funcs = default_funcs(p.numinputs)
  ev_count = 0
  for _ = 1:nreps
    # evolve a chromosome that outputs goal g.
    res = mut_evolve_repeat( 10, p, [g], funcs, maxsteps)
    if res == nothing
      error("run_random_neutral_walk: mut_evolve_repeat() failed.  Try increasing maxsteps.")
    end
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = res
    #(c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve_repeat( 10, p, [g], funcs, maxsteps)
    @assert output_values(c) == g
    ev_count += random_neutral_walk( c, steps, maxsteps, skip_steps=skip_steps, evo_count_only=true )
  end
  ev_count/nreps
end

function run_random_neutral_walk( c_list::Vector{Chromosome}, steps::Int64, maxsteps::Int64, nreps::Int64;
    csvfile::String="" )
  df = DataFrame()
  df.complexity = Float64[]
  df.avg_complexity = Float64[]
  df.all_count = Float64[]
  df.robust_count = Float64[]
  df.evolvable_count = Float64[]
  df.evolvability = Float64[]
  df.circuit_dist = Float64[]
  println("size(df): ",size(df))
  p = c_list[1].params
  for c in c_list
    @assert c.params == p
  end
  row_avg_list = pmap( c->avg_random_neutral_walks( c, steps, maxsteps, nreps ), c_list )
  #row_avg_list = map( c->avg_random_neutral_walks( c, steps, maxsteps, nreps ), c_list )
  for row_avg in row_avg_list
    push!(df,row_avg)
  end
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  println("# host: ",hostname," with ",nprocs()-1,"  processes: " )
  println("# date: ",Dates.now())
  println("# funcs: ", funcs)
  println("# steps: ",steps)
  println("# maxsteps: ",maxsteps)
  println("# nreps: ",nreps)
  open( csvfile, "w" ) do f
    println(f,"# date : ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    println(f,"# funcs: ", funcs)
    print_parameters(f,p,comment=true),
    println(f,"# steps: ",steps)
    println(f,"# maxsteps: ",maxsteps)
    println(f,"# nreps: ",nreps)        
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

# Helper function for the above run_random_neutral_walk()
function avg_random_neutral_walks( c::Chromosome, steps::Int64, maxsteps::Int64, nreps::Int64 )
    df_row = Tuple[]
    for i = 1:nreps
      push!(df_row, random_neutral_walk( c, steps, maxsteps, evo_count_only=false ) )
      #println("df_row: ",df_row[end])
    end
    df_row_avg = [ mean(df_row[k][j] for k = 1:length(df_row)) for j = 1:length(df_row[end]) ]
    #println("df_row_avg: ",df_row_avg)
  df_row_avg
end 
      
# Returns an evolvable_count or a dataframe row depending on whether evo_count_only is true or false
function random_neutral_walk( c::Chromosome, steps::Int64, maxsteps::Int64; evo_count_only::Bool=false )
  orig_complexity =  complexity5(c)
  c_orig = deepcopy(c)
  all_outputs_list = Goal[]
  all_unique_list = Goal[]
  avg_complexity = zeros(Float64,steps)   # mean complexity of c over steps
  all_count = zeros(Int64,steps)
  robust_count = zeros(Int64,steps)
  evolvable_count = zeros(Int64,steps)
  evolvability = zeros(Float64,steps)
  circuit_dist = zeros(Float64,steps)
  evolvable_count = zeros(Int64,steps)
  evolvability = zeros(Float64,steps)
  funcs = default_funcs(c.params.numinputs)
  goal = output_values(c)
  total_steps = 0  # Number of mutation steps 
  #println("goal: ",goal,"  complexity: ",complexity5(c))
  for i = 1:steps
    j = 1
    while j <= maxsteps   # terminated by a break statement
      sav_c = deepcopy(c)
      (c,active) = mutate_chromosome!( c, funcs )
      outputs = output_values(c)
      if outputs == goal
        #println("successful step for i= ",i,"  outputs: ",outputs)
        break
      end
      c = sav_c
      j += 1
    end
    if step == maxsteps
      error("Failed to find neutral mutation")
    end
    total_steps += j
    avg_complexity[i] = complexity5(c)
    all_outputs = mutate_all( c, funcs )
    robustness = length(filter(x->x==goal,all_outputs))
    robust_count[i] =  i == 1 ? robustness : robust_count[i-1] + robustness
    all_count[i] =  i == 1 ? length(all_outputs) : all_count[i-1] + length(all_outputs)
    all_unique_list = unique( vcat( all_unique_list, unique( all_outputs)))
    evolvable_count[i] = length(all_unique_list)
    evolvability[i] = evolvable_count[i]/all_count[i]
    circuit_dist[i] = circuit_distance( c, c_orig )
    #println("k: ",k,"  all_count[i]: ",all_count[i],"  robust_count[i]: ",robust_count[i],"  evolvable_count[i]: ",evolvable_count[i] )
  end
  println("total steps: ",total_steps)
  if evo_count_only
    evolvable_count[i]
  else
    ( orig_complexity, avg_complexity[steps], all_count[steps], robust_count[steps], evolvable_count[steps], 
      evolvability[steps], circuit_dist[steps] )
  end
end

function binned_circuit_complexities( p::Parameters, n_circuits::Int64, max_bin_size::Int64, 
    complexity_lower_bounds::Vector{Float64} )
  funcs = default_funcs(p.numinputs)
  complexity_bins = [ Chromosome[] for k = 1:length(complexity_lower_bounds) ]
  for i = 1:n_circuits 
    c = random_chromosome( p, funcs )
    cmplx = complexity5(c)
    cmplx = cmplx >= 0.0 ? cmplx : 0.0
    k = 1
    while k < length(complexity_lower_bounds) && cmplx >= complexity_lower_bounds[k+1]
      k += 1 
    end
    #println("cmplx: ",cmplx,"  k: ",k)
    if length(complexity_bins[k]) < max_bin_size
      push!(complexity_bins[k], c )
    end
  end
  complexity_bins
end
