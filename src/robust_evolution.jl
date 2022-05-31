# First written on 5/11/21.  Finally debugged on 5/13/21.
# Given a target phenotype, employs neutral (phenotype doesn't change) evolution to discover increeasingly robust genotypes that map to the target phenotype.
#using HypothesisTests
using Random


# Returns result_list which is a list of (circuit,rbst) pairs where rbst is the robustness() of the circuit.
# The first element of the list has maximum robustness.
function robust_evolution_parallel( p::Parameters, phlist::GoalList, nreps::Int64, max_steps::Integer, max_mutates::Int64, max_tries::Int64, funcs::Vector{Func};
      use_neutral::Bool=true, use_lincircuit::Bool=false, csvfile::String="" )
  if use_neutral
    #result_list = pmap( ph->robust_evolution_neutral( p, ph, nreps, max_steps, max_tries, funcs, use_lincircuit=use_lincircuit ), phlist )
    result_list = map( ph->robust_evolution_neutral( p, ph, nreps, max_steps, max_tries, funcs, use_lincircuit=use_lincircuit ), phlist )
  else
    #result_list = pmap( ph->robust_evolution( p, ph, nreps, max_steps, max_mutates, max_tries, funcs, use_lincircuit=use_lincircuit ), phlist )
    result_list = map( ph->robust_evolution( p, ph, nreps, max_steps, max_mutates, max_tries, funcs, use_lincircuit=use_lincircuit ), phlist )
  end
  df = DataFrame()
  df.goal = String[]
  df.avg_start_rbst = Float64[]
  df.best_rbst = Float64[]
  for i = 1:length(phlist)
    goal = @sprintf("0x%04x",phlist[i][1])
    startrbst = result_list[i][2]
    bestrbst = result_list[i][3]
    println("i: ",i," ",goal," startrbst: ",startrbst,"  bestrbst: ",bestrbst)
    push!(df,(goal,startrbst,bestrbst))
  end
  if length(csvfile) > 0
    write_df_to_csv(df,p,funcs,csvfile,nreps=nreps,max_steps=max_steps,max_mutates=max_mutates)
  end
  df
end

# This version evolves many genotypes that map to the given phenotype ph and thus is not an example of doing a random phenotype-preserving walk starting at a given genotype
# returns a (circuit,rbst) pair where rbst is the robustness() of the circuit.
function robust_evolution( p::Parameters, ph::Goal, nreps::Int64, max_steps::Integer, max_mutates::Int64, max_tries::Int64, funcs::Vector{Func}; use_lincircuit::Bool=false )
  circ_rbst_list = Tuple{Circuit,Float64}[]
  sum_start_rbst = 0.0
  for rep = 1:nreps
    #println("rep: ",rep)
    tries = 0
    done = false
    nc = c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
    while !done && tries < max_tries
      (nc,steps) = neutral_evolution( c, ph, max_steps )
      tries += 1
      if steps < max_steps
        done = true
      end
      c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
    end
    if !done 
      error("neutral evolution failed in function robust_evolution() ")
    end
    @assert ph == output_values(nc)
    rbst = robustness(nc,funcs)
    sum_start_rbst += rbst
    #print("step: ",0,"  ph: ",ph,"  rbst: ",rbst,"   ")
    #print_circuit(c,funcs)
    for step = 1:max_mutates
      #(new_c,active) = mutate_chromosome!( deepcopy(nc), funcs )
      if use_lincircuit
        new_c = mutate_circuit!(deepcopy(c),funcs)
      else
        (new_c,active) = mutate_chromosome!(deepcopy(nc),funcs)
      end
      new_rbst = robustness( new_c, funcs )
      new_ov = output_values( new_c )
      #println("step: ",step,"  new_ov: ",new_ov,"  new_rbst: ",new_rbst)
      if new_ov != ph || new_rbst < rbst
        continue
      end
      #println("step: ",step,"  new ov: ",output_values(new_c),"  new_rbst: ",new_rbst)
      c = new_c
      @assert output_values(c) == ph 
      rbst = new_rbst
      #print("step: ",step,"  ov: ",output_values(c),"  rbst: ",rbst,"   ")
      #print_circuit(c,funcs)
    end
    push!(circ_rbst_list, (c,rbst))
  end
  average_starting_robustness = sum_start_rbst/nreps
  sort!(circ_rbst_list,by=x->x[2])  # sort on robustness
  (circ,best_rbst) =  circ_rbst_list[end]        # return the pair with the largest robustness
  (circ,average_starting_robustness,best_rbst)
end

# This version does a random phenotype preseriving walk starting with a single circuit that maps to the given phenotype ph 
# returns a (circuit,rbst) pair where rbst is the robustness() of the circuit. 
function robust_evolution_neutral( p::Parameters, ph::Goal, nreps::Int64, max_steps::Integer, max_tries::Int64, funcs::Vector{Func}; use_lincircuit::Bool=false )
  tries = 0
  done = false
  nc = c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
  while !done && tries < max_tries
    #println("ph: ",ph,"  c: ",c)
    (nc,steps) = neutral_evolution( c, ph, max_steps )
    tries += 1
    if steps < max_steps
      done = true
    end
    c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
  end
  if !done 
    error("neutral evolution failed in function robust_evolution() ")
  end
  @assert ph == output_values(nc)
  rbst = starting_robustness = robustness(c,funcs)
  #print("step: ",0,"  ph: ",ph,"  rbst: ",rbst,"   ")
  print_circuit(c,funcs)
  for rep = 1:nreps
    #(new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    if typeof(c)==LinCircuit
      new_c = mutate_circuit!(deepcopy(c),funcs)
    else
      (new_c,active) = mutate_chromosome!(deepcopy(c),funcs)
    end
    new_rbst = robustness( new_c, funcs )
    new_ov = output_values( new_c )
    #println("rep: ",rep,"  new ov: ",output_values(new_c),"  new_rbst: ",new_rbst)
    if new_ov != ph || new_rbst < rbst
      continue
    end
    #println("rep: ",rep,"  new ov: ",output_values(new_c),"  new_rbst: ",new_rbst)
    c = new_c
    @assert output_values(c) == ph
    rbst = new_rbst
    print("step: ",step,"  ov: ",output_values(c),"  rbst: ",rbst,"   ")
    #print_circuit(c,funcs)
  end
  (c,starting_robustness,rbst)
end

# Not used
function robust_evolution_csv( c::Circuit, max_steps::Integer, nreps::Int64, funcs::Vector{Func}, csvfile::String )
  df = DataFrame()
  df.numinputs = Int64[]
  df.numints = Int64[]
  df.numlb = Int64[]
  df.max_steps = Int64[]
  df.nreps = Int64[]
  df.goal = Vector{MyInt}[]
  df.start_rbst = Float64[]
  df.final_rbst = Float64[]
  g = output_values(c,funcs)
  start_rbst = robustness(c,funcs)
  rlist = robust_evolution_parallel(c,max_steps,nreps,funcs)
  push!(df,(c.params.numinputs,c.params.numinteriors,c.params.numlevelsback,max_steps,nreps,g,start_rbst,rlist[1][2]))
  write_df_to_csv( df, p, funcs, csvfile, max_steps=max_steps, nreps=nreps )
  df
end
  
# g is the target phenotype.  Not used
function robust_evolution_pop( c::Circuit, popsize::Integer, max_steps::Integer, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs=default_funcs(c.params)
  end
  g = output_values( c )
  print("g: ",g,"  popsize: ",popsize,"  max_steps: ",max_steps," length(funcs): ",length(funcs),"    ")
  print_circuit(c,funcs)
  rbst = robustness(c,funcs)
  ov = g
  print("ov: ", ov, "  rbst: ", rbst, "   " )
  print_circuit(c,funcs)
  for step = 1:max_steps 
    pop = Tuple{Chromosome,Float64}[]
    for i = 1:popsize
      #(new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      if typeof(c)==LinCircuit
        new_c = mutate_circuit!(deepcopy(c),funcs)
      else
        (new_c,active) = mutate_chromosome!(deepcopy(c),funcs)
      end
      new_ov = output_values( new_c )
      new_rbst = robustness( new_c, funcs )
      if new_ov == ov && new_rbst >= rbst
        println("step: ",step,"  new ov: ",output_values(new_c),"  new_rbst: ",new_rbst)
        push!(pop,( deepcopy(new_c), new_rbst ))
      end
    end
    #print_pop(pop)
    if length(pop) > 1
      sort!(pop,lt=my_isless,rev=true)
    end
    if length(pop) > 0
      println("new c:  ")
      print_pop(pop)
      c = pop[1][1]
      rbst = pop[1][2]
    end
  end
  println("robust_evolution returned with phenotype: ", ov, " with robustness: ", rbst )
  @assert ov == g
  (c,rbst)
end
      
function print_pop( pop::Vector{Tuple{Chromosome,Float64}} )
  for i = 1:length(pop)
    #print( "ov: ",@sprintf("0x%04x",pop[i][2][1]),"  rbst: ",pop[i][2], "    " )
    print( "ov: ",output_values(pop[i][1]),"  rbst: ",pop[i][2], "    " )
    print_circuit(pop[i][1],funcs)
  end
end

function print_rlist( rlist )
  for r in rlist
    print("r[2]: ",r[2]," ov: ",output_values(r[1],funcs),"   ")
    print_circuit(r[1],funcs)
  end
end

# Compare x and y on the second component
function my_isless( x::Tuple, y::Tuple )
  return x[2] < y[2]
end

#=
function mutate_circuit!( c::Circuit, funcs::Vector{Func} )
  if typeof(c) == Chromosome
    mutate_chromosome!( c, funcs )
  elseif typeof(c) == LinCircuit
    mutate_circuit!( c, funcs )
  end
end
=#    

function robustness( c::Circuit, funcs::Vector{Func} )
  #print("robustness: c:  ")
  #print_circuit(c,funcs)
  c_output = output_values(c,funcs)
  outputs = mutate_all( c, funcs, output_outputs=true )
  #println("outputs[1]: ",outputs[1])
  robust_outputs = filter( x->x==c_output, outputs )
  return length(robust_outputs)/length(outputs)
end     
