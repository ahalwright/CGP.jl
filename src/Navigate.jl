# Attempts to compute a fitness non-decreasing path from the circuit represented by circ_int to a genotype of the phenotype phdest.
# Sometimes works, but often goes into a long loop.  
# However, our exact evolvability results from the evolvability paper show that there is a mutation between any pair of phenotypes
#    for the 3-input 8-gates case, which proves perfect navigability.
using DataStructures

# Creates a dataframe with frequency and K complexity properties of runs of function epochal_evolution_fitness().  
# If one_success == true, terminate iteration over circuits when a success of epochal_evolution_fitness is found.
# If csvfile is a legitimate filename, writes this dataframe to this file.
# function dataframe_row() is pmap called ngoals times.  dataframe_row() calls function navigate().
function run_navigate_epochal( P::Parameters, funcs::Vector{Func}, fitness::Vector{Float64}, ngoals::Int64, max_steps::Int64, num_circuits::Int64; 
      use_pmap::Bool=false, one_success::Bool=false, csvfile::String="" )
  rdict = redundancy_dict(P,funcs)
  kdict = kolmogorov_complexity_dict(P,funcs)
  rdf = DataFrame( :numinputs=>Int64[], :numgates=>Int64[], :srcfit=>Float64[], :srcK=>Int64[], :srcfreq=>Int64[], 
      :destfit=>Float64[], :destK=>Int64[], :destfreq=>Int64[], :mean_steps=>Float64[], :num_circuits=>Int64[], :ncircuits=>Int64[], :failures=>Int64[]  )
  if use_pmap
    row_list = pmap( _->dataframe_row( P, funcs, fitness, kdict, rdict, num_circuits, max_steps, one_success=one_success ), 1:ngoals )
  else
    row_list = map( _->dataframe_row( P, funcs, fitness, kdict, rdict, num_circuits, max_steps, one_success=one_success ), 1:ngoals )
  end
  for row in row_list
    push!( rdf, row )
  end
  count_failures = length(findall( x->x==Float64(max_steps),rdf.mean_steps))
  println("count_failures: ",count_failures)
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,P,comment=true)
      println(f,"# funcs: ", funcs )
      println(f,"# max_steps: ",max_steps)
      println(f,"# num_circuits: ",num_circuits)
      println(f,"# ngoals: ",ngoals)
      println(f,"# one_success: ", one_success ? "true" : "false" )
      println(f,"# count failures: ",count_failures)
      CSV.write( f, rdf, append=true, writeheader=true )
    end
  end
  rdf
end

# Chooses a random destgoal and then chooses a random srcgoal with smaller fitness than destgoal.
# Returns a tuple which becomes a row of the dataframe defined in function run_navigate_epoch().
function dataframe_row( P::Parameters, funcs::Vector{Func}, fitness::Vector{Float64}, kdict::Dict{MyInt,Int64}, rdict::Dict{MyInt,Int64},
        num_circuits::Int64,  max_steps::Int64; one_success::Bool=false )
  #( msteps, failures ) = navigate( P, funcs, fitness, srcgoal, destgoal, max_steps, num_circuits=num_circuits, use_pmap=false )
  fitfunct( g::Goal ) = fitness[ g[1]+1 ]
  quantile_min = 0.1
  destgoal = randgoal(P)
  while fitfunct( destgoal ) < quantile(fitness,quantile_min)
    destgoal = randgoal(P)
  end
  srcgoal = randgoal(P)
  while fitfunct( srcgoal ) > fitfunct( destgoal )  # Choose srcgoal to have fitness less than or equal to destgoal
    srcgoal = randgoal(P)
  end
  println("srcgoal: ",srcgoal,"  destgoal: ",destgoal)
  ( msteps, failures, ncircuits ) = navigate( P, funcs, fitness, srcgoal, destgoal, max_steps, num_circuits=num_circuits, one_success=one_success, use_pmap=false )
  println("(msteps,failures): ",(msteps,failures))
  return ( P.numinputs, P.numinteriors, fitfunct(srcgoal), kdict[srcgoal[1]], rdict[srcgoal[1]], fitfunct(destgoal), kdict[destgoal[1]], rdict[destgoal[1]], msteps, num_circuits, ncircuits, failures  )
end

# fitness is a given vector over phenotypes with values in the interval [0,1].
# goal1 and goal2 are phenotypes which are Vectors of MyInts (one MyInt per output).
# If fitness of goal1 is greater than fitness of goal2, then goal1 and goal2 are switched. Then srcgoal is goal1, destgoal is goal2.
# Then srcgaol is the lower fitness goal and destgoal is the higher fitness goal.
# Circuits that map to srcgoal are found by calling neutral_evolution() (which is really epochal_evolution).
# This function calls function epochal_evolution_fitness() on num_circuits circuits evolved to map to srcgoal.
# Returns the pair of the mean number of steps (including for unsuccessful runs) and the number of failed runs of epochal_evolution_fitness().
# Note: If navigate() is called in a pmap(), use_pmap must be false because nesting pmaps does not work.
function navigate( P::Parameters, funcs::Vector{Func}, fitness::Vector{Float64}, goal1::Goal, goal2::Goal, max_steps::Int64;
      num_circuits::Int64, one_success::Bool=false, use_pmap::Bool=false )::Union{Tuple{Float64, Int64, Int64},Nothing}
  println("function navigate.  one_success: ", one_success )
  ncircuits = num_circuits
  max_tries = 10
  srcgoal = goal1
  destgoal = goal2
  if fitness[ srcgoal[1]+1 ] > fitness[destgoal[1]+1]
    tmpgoal = deepcopy(destgoal)
    destgoal = deepcopy(srcgoal)
    srcgoal = tmpgoal
  end
  println("srcgoal: ",srcgoal,"  fitness(srcgoal): ",fitness[srcgoal[1]+1], "  destgoal: ",destgoal,"  fitness(destgoal): ",fitness[destgoal[1]+1] )
  steps_list = Int64[]
  for i = 1:num_circuits
    rch = random_chromosome( P, funcs )
    (cch, step) = neutral_evolution( rch, funcs, srcgoal, max_steps, print_steps=false )
    ttry = 0
    while ttry < max_tries &&  step==max_steps
      (cch, step) = neutral_evolution( rch, funcs, srcgoal, max_steps, print_steps=false )
      ttry += 1
    end
    if step == max_steps
      println("neutral evolution failed in function navigate().  i: ",i)
      return nothing
    end
    (ach, steps, src_fit, dest_fit ) = epochal_evolution_fitness( cch, funcs, destgoal, fitness, max_steps, print_steps=false )
    println("i: ",i,"  epochal_evol_fitness: steps: ",steps)
    push!(steps_list,steps)
    if one_success && steps < max_steps
      ncircuits = i
      println("one_success break.  ncircuits: ",ncircuits)
      break   # exit for i = 1:num_circuits loop
    end
  end  
  failures = length( findall(x->x==max_steps,steps_list))
  return ( mean( steps_list ), failures, ncircuits )
end


#= obsolete
function run_navigate_both( P::Parameters, funcs::Vector{Func}, fitness::Vector{Float64}, ngoals::Int64, max_steps::Int64, max_reps::Int64;
    num_circuits::Int64=10, use_pmap::Bool=false, suffix::String="X" )
  df_name = "count_outputs_ch_$(length(funcs))funcs_$(P.numinputs)inputs_$(P.numinteriors)gates_$(P.numlevels back)lb_$(suffix).csv"
  df = read_dataframe( "../data/counts/$(df_name)")
  fitfunct( g::Goal ) = fitness[ g[1]+1 ]   # local function
  rdict = redundancy_dict(P,funcs)
  kdict = kolmogorov_complexity_dict(P,funcs)
  quantile_min = 0.1 
  rdf = DataFrame( :numinputs=>Int64[], :numgates=>Int64[], :srcfit=>Float64[], :srcK=>Int64[], :srcfreq=>Int64[], 
      :destfit=>Float64[], :destK=>Int64[], :destfreq=>Int64[], :failures=>Tuple{Int64,Int64}[], :mean_steps=>Float64[], 
      :bfs_failures=>Int64[], :bfs_sum_reps=>Int64[] )
  for i = 1:ngoals
    destgoal = randgoal(P)
    while fitfunct( destgoal ) < quantile(fitness,quantile_min)
      destgoal = randgoal(P)
    end
    srcgoal = randgoal(P)
    while fitfunct( srcgoal ) > fitfunct( destgoal )  # Choose srcgoal to have fitness less than or equal to destgoal
      srcgoal = randgoal(P)
    end
    src_circuits = string_to_expression( df.circuits_list[srcgoal[1]+1] )
    ( msteps, failures ) = navigate( P, funcs, fitness, srcgoal, destgoal, max_steps, df, num_circuits=num_circuits, use_pmap=use_pmap )
    println("(msteps,failures): ",(msteps,failures))
    sum_reps = 0
    bfs_failures = 0
    for i = 1:num_circuits
      ( reps, reverse_fits ) = navigate_bfs( P, funcs, src_circuits[i], fitness, destgoal, max_reps=max_reps )
      sum_reps += reps
      bfs_failures = reps==max_reps ? bfs_failures+1 : bfs_failures
    end
    println("bfs_failures: ",bfs_failures)
    println("( reps, reverse_fits ): ",( reps, reverse_fits ))  
    push!( rdf, 
      ( P.numinputs, P.numinteriors, fitfunct(srcgoal), kdict[srcgoal[1]], rdict[srcgoal[1]], fitfunct(destgoal), kdict[destgoal[1]], rdict[destgoal[1]], 
          failures, msteps, bfs_failures, sum_reps ) )
  end
  rdf
end
=#

#= Use breadth-first search to find a fitness nondecreasing path
function find_nondecreasing_BFS( ( c::Circuit, funcs::Vector{Func}, g::Goal, fitfunct::Union{Vector{Float64},Function}, max_steps::Integer;
   print_steps::Bool=false )::Tuple{Union{Circuit,Nothing},Int64,Float64,Float64,Float64}
end
=#

# Example run:  run_navigate_bfs( P3, funcs, fitness, 2, 10_000, print_interval=1000 )
function run_navigate_bfs( P::Parameters, funcs::Vector{Func}, fitness::Vector{Float64}, ngoalpairs::Int64, max_reps::Int64, max_steps::Int64=100_000; print_interval::Int64=1000, 
    csvfile::String="" )
  function nbfs_helper( df::DataFrame )
    goal1 = randgoal(P)
    goal2 = randgoal(P)
    if fitness[goal1[1]+1] <= fitness[goal2[1]+1]
      srcgoal = goal1
      destgoal = goal2
    else
      srcgoal = goal2
      destgoal = goal1
    end
    ( ch, steps ) = neutral_evolution( random_chromosome(P,funcs), funcs, srcgoal, max_steps )
    ( reps, path ) = navigate_bfs( P, funcs, ch, fitness, destgoal, max_reps=max_reps, print_interval=print_interval )
    println("reps: ",reps,"  path: ",path)
    path_length = path==nothing ? 0 : length(path)
    ( srcgoal, destgoal, max_reps, reps, path_length ) 
  end
  df = DataFrame( :srcgoal=>Goal[], :destgoal=>Goal[], :max_reps=>Int64[], :reps=>Int64[], :path_length=>Int64[] )
  row_list = pmap( _->nbfs_helper(df), 1:ngoalpairs )
  for row in row_list
    push!(df, row )
  end
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,P,comment=true)
      println(f,"# funcs: ", funcs )
      println(f,"# max_steps: ",max_steps)
      println(f,"# max_reps: ",max_reps)
      println(f,"# ngoalpairs: ",ngoalpairs)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function navigate_bfs( P::Parameters, funcs::Vector{Func}, ch::Chromosome, fitness::Vector{Float64}, phdest::Vector{MyInt}; max_reps::Int64=10, print_interval::Int64=1000 )
  #src_ph = output_values(int_to_chromosome( circ_int, P, funcs ))
  circ_int = chromosome_to_int( ch, funcs )
  src_ph = output_values(ch)
  src_fit = fitness[ src_ph[1] + 1 ]
  dest_fit = fitness[ phdest[1] + 1 ]      
  println("(src_fit,dest_fit): ",(src_fit,dest_fit))
  @assert src_fit <= dest_fit
  prev_circuit_dict = Dict{Int128,Tuple{Int128,Int64}}()
  #explored = Set{Int128}( [] )
  queue = Queue{Int128}()    # A circuit_int and its predecessor
  prev_circuit_dict[circ_int] = (circ_int,0)
  enqueue!( queue, circ_int )
  #println("queue: ", queue)
  currrent_fit = src_fit
  reps = 0
  while length(queue) > 0 && reps < max_reps
    reps += 1
    if reps % print_interval == 0
      println("while reps: ",reps)
    end
    current_int = dequeue!( queue )
    current_circ = int_to_chromosome( current_int, P, funcs )
    current_ph = output_values(current_circ)
    current_fit = fitness[ current_ph[1] + 1 ]
    (prev_key,current_plength) = prev_circuit_dict[ current_int ]
    #println("current_int: ",current_int,"  current_plength: ",current_plength,"  current_ph: ",current_ph,"  current_fit: ",current_fit,"  prev_key: ",prev_key)  
    #push!( explored, current_int )
    mut_circs = mutate_all( current_circ, funcs, output_outputs=false, output_circuits=true )
    for mcirc in mut_circs
      mph = output_values(mcirc)
      mfit = fitness[ mph[1] + 1 ]      
      mcirc_int = chromosome_to_int( mcirc, funcs )
      #println("mph: ",mph,"  mfit: ",mfit,"  mcirc_int: ",mcirc_int)
      if mph == phdest   # Goal found, produce result
        println("found phdest")
        reverse_path = Int128[mcirc_int,current_int]
        done = false
        cnt = 0
        cur_int = current_int
        while haskey( prev_circuit_dict, cur_int ) && !done && cnt < 10
          ( prev_int, prev_plength ) = prev_circuit_dict[ cur_int ]
          #println("( prev_int, prev_plength ): ",( prev_int, prev_plength ))
          push!( reverse_path, prev_int )
          if prev_int == current_int || prev_int == circ_int 
            done = true
          end
          cur_int = prev_int
          cnt += 1; # println("cnt increased")
        end
        reverse_fits = map( x->(x,fitness[ output_values( int_to_chromosome( x, P, funcs ))[1]+1]), reverse_path )
        #println("reverse_fits: ",reverse_fits)
        return (reps,reverse_fits)
      end
      if (mfit <= dest_fit && mfit >= current_fit)
        if haskey( prev_circuit_dict, mcirc_int)
          ( save_prev_key, saved_plength ) = prev_circuit_dict[ mcirc_int ]
          #println("haskey ( save_prev_key, saved_plength ): ",( save_prev_key, saved_plength ))
          if current_plength + 1 < saved_plength
            prev_circuit_dict[ mcirc_int ] = (current_circuit, current_plength + 1)
            println("revise dict: (current_circuit, current_plength + 1): ",(current_circuit, current_plength + 1))
          end
        else
          prev_circuit_dict[ mcirc_int ] = (current_int, current_plength + 1)
          enqueue!( queue, mcirc_int )
          #if length(queue) % 100 == 0
          #  println("length(queue): ", length(queue))
          #end
        end
      end
    end  # for loop
  end # while loop
  println("function navigate_bfs() finished with an empty queue.")
  return (reps,nothing)
end
  
#=
function navigate_bfs( P::Parameters, funcs::Vector{Func}, circ_int::Int128, fitness::Vector{Float64}, phdest::Vector{MyInt} )
  @assert length(fitness) == 2^2^P.numinputs
  circ = int_to_chromosome( circ_int, P, funcs )
  phsrc = output_values( circ )
  src_fit = fitness[ phsrc[1]+1 ]
  println("phsrc: ",phsrc,"  src fit: ",src_fit)
  println("phdest: ",phdest,"  dest fit: ",fitness[ phdest[1]+1 ])
  @assert src_fit <= fitness[ phdest[1]+1 ]
  #path = Tuple{Int128,Float64}[(circ_int,fit)]
  queue = Queue{Tuple{Int128,Int128}}()    # A circuit_int and its predecessor
  #explored = Set{Vector{Tuple{Int128,Float64}}}( [] )
  explored = Set{Tuple{Int128,Int128}}( [] )
  current_pair = (circ_int,circ_int)
  enqueue!( queue, current_pair )
  while length(queue) > 0
    #println("queue: ", queue )
    current_pair = dequeue!( queue )
    push!( explored, current_pair )
    current_circ = int_to_chromosome( current_pair[1], P, funcs )
    prev_circ = int_to_chromosome( current_pair[2], P, funcs )
    print("current_circ:  ")
    print_circuit(current_circ )
    print("prev_circ:  ")
    print_circuit(prev_circ)
    #println("prev fit: ")
    mut_circs = mutate_all( current_circ, funcs, output_outputs=false, output_circuits=true )
    for mcirc in mut_circs
      mcirc_int = chromosome_to_int( mcirc, funcs )
      prev_int = 
      if !(mcirc_int in explored)
        mph = output_values( mcirc )
        fit = fitness[ mph[1]+1 ]
        #print("mcirc:  fit: ",fit,"  mph: ",mph,"   ")
        #print_circuit(mcirc)
        if fitness[phsrc[1]+1] <= fit && fit <= fitness[phdest[1]+1] && fit >= current_path[end][2]
          mcirc_int = chromosome_to_int( mcirc, funcs )
          mpath = push!(path,(mcirc_int,fit))
          println("mph: ",mph,"  mcirc_int: ",mcirc_int,"  phdest: ",phdest)
          if mph == phdest
            return mpath
          end
          #if !(mpath in explored) 
          enqueue!( queue, mpath )
          push!( explored, mcirc_int )
          #end
        end
      end
    end
  end
end
=#
