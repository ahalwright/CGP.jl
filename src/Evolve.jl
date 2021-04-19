export components_matched, goals_matched_exact, goals_matched_hamming, next_chromosome!, mut_evolve, mut_evolve_increase_numints
export mut_reduce_numactive
export random_neutral_walk, match_score, findmaxall, findminall, findmaxall_1, findminall_1, findmaxrand, findminrand
export evolve_function, mut_evolve_repeat, circuit_evolve, run_circuit_evolve
export neutral_evolution, geno_circuits, geno_properties, geno_list_properties

# 5/21: To test:
# ni = 2; nc = 4
# julia> g = randgoallist(2,ni,nc)  # For this example g[1] is interpreted as output, g[2] is goal
# julia> pmin=match_score(g[1],g[2],ni)
# Not efficient if length(output) is greater than about 6
# Returns score of the best-scoring match between components of output and components of goal.
# A match between the components of outputs and the components of goal is permutation p of 1:nc
#   where ouput[i] is matched to goal[p[i]].
# There are two options depending on the keyword parameter avgfitness.
# match_score is the number of exact component matches plus the score of the best partial component match.
#
# There are two options depending on the setting of avgfitness.
# The score of of a partial match between component x and component y is 1.0 - hamming(x,y)/2^numinputs.
# Example:  let output = [0x0e, 0x0d, 0x00, 0x04], goal = [0x06, 0x04, 0x0a, 0x01]
# There is an exact match between ouput[4] and goal[2].
# One of the best partial matches is between output[1] and goal[1] where hamming(output[1],goal[1]) = 1.
# This gives a partial match score of 0.75.
# Thus, match_score(output,goal,ni) returns 1.75
function match_score( output::Vector{MyInt}, goal::Vector{MyInt},numinputs::Int64;
      avgfitness::Bool=false)
  #println("match_score: avgfitness: ",avgfitness)
  nc = length(output)  # number of components
  #println("output: ",output,"  goal: ",goal)
  @assert nc == length(goal)
  H = [ hamming_distance(output[i],goal[j],numinputs) for i = 1:nc, j=1:nc]  
  #println("H: ",H)
  P = permutations(collect(1:nc))
  if avgfitness
    mxscores = [ (sum( 1.0-H[i,p[i]] for i = 1:nc ),[p[i] for i = 1:nc]) for p in P ]
    #println("mxscores: ",mxscores)
    findmaxall(mxscores)[1][1]
  else # number of exact matches plus Hamming score for best partial match
    mxscores = [ maxscore([H[i,p[i]] for i = 1:nc]) for p in P ]
    findmaxall(mxscores)[1]
  end
end
#=
function match_score( output::Vector{MyInt}, goal::Vector{MyInt},numinputs::Int64)
  nc = length(output)  # number of components
  println("match_score: output: ",output,"  goal: ",goal)
  @assert nc == length(goal)
  H = [ hamming_distance(output[i],goal[j],numinputs) for i = 1:nc, j=1:nc]  
  println("H: ",H)
  P = permutations(collect(1:nc))
  mxscores = [ maxscore([H[i,p[i]] for i = 1:nc]) for p in P ]
  println("mxscores: ",mxscores)
  findmaxall(mxscores)[1]
end
=#

# return the number of zeros plus one minus the smallest non-zero value
# Thus, if V contains 2 zeros and the smallest non-zero value is 0.2, returns 2.8.
# If the components of V are Hamming distances between components x and y, then 
#   maxscore is a rating of the quality of match between x and y
function maxscore( V::Vector{Float64} )
  minind = 1
  count_zeros = 0
  for i = 1:length(V)
    if V[i] > 0 
      if (V[minind] <= 0 || V[i] < V[minind])
        minind = i
      end
    else
      count_zeros += 1
    end
  end
  result = count_zeros < length(V) ? (count_zeros + 1.0 - V[minind]) : Float64(count_zeros)
  @assert result <= length(V)
  result
end

# Returns a list of triples (comp,i,j) where comp==output[i]==goal[j] and for different triples
#    (comp1,i1,j1) and (comp2,i2,j2) we must have i1!=i2 and j1!=j2 (but we may have comp1==comp2).
# So the length of the result is the degree of matching between outputs and goal
# Based on a merge sort algorithm
# EXample:  let output=[0x7,0xd], let goal = [0x1,0x7]
# Result (0x07, 1, 2)  which says that output[1] == goal[2] == 0x07
function components_matched( output::Vector{MyInt}, goal::Vector{MyInt} )
  result = Tuple{MyInt,Int64,Int64}[]
  len = length(output)
  output_pairs = [(output[i],i) for i = 1:len]
  goal_pairs = [(goal[i],i) for i = 1:len]
  spoutput = sortperm(output_pairs)   # permutation that sorts the output components
  spgoal = sortperm(goal_pairs)         # permutation that sorts the goal components
  #println("spoutput: ",spoutput)  
  #println("spgoal: ",spgoal)
  #println("soutput: ",output[spoutput])
  #println("sgoal: ",goal[spgoal])
  i = 1
  j = 1
  while i <= len && j <= len
    if output[spoutput[i]][1] == goal[spgoal[j]][1]
      push!(result,(output_pairs[spoutput[i]][1],output_pairs[spoutput[i]][2],goal_pairs[spgoal[j]][2]))
      i += 1
      j += 1
    elseif output[spoutput[i]][1] < goal[spgoal[j]][1]
      #println()
      i += 1
    elseif output[spoutput[i]][1] > goal[spgoal[j]][1]
      #println()
      j += 1
    end
  end
  result
end

# returns a triple (num_components_matched, goals_matched_list, components_matched_list)
#  where 
#   num_components_matched is the maximum number of components of output which match a goal
#   goals_matched_list is the list of goals that achieve this maximum
#   components_matched_list is the list of components matched for each of the goals that achieves the maximum
#       (these are the outputs of the components_matched() function for each goal that achieves the maximum)
# Example:  
#  output = [0x7,0xd]
#  goallist = [[0x1,0x2],[0x1,0x7],[0x7,0x3]]
#  goals_matched(output,goallist) = (1, [2, 3], Array{Tuple{UInt8,Int64,Int64},1}[[(0x07, 1, 2)], [(0x07, 1, 1)]])
#    which says that at most 1 component of output matches with a goal, and this is achieved for goal[2] and goal[3].
#    See the comments for components_matched() for the interpreation of the third component of the result triple.
# Note that fitness does not take robustness into account
function goals_matched_exact( output::Vector{MyInt}, goallist::GoalList, numinputs::Int64;
      avgfitness::Bool=false, perm_heuristic::Bool=false )   # Optional arguments not used in this function
  #println("exact")
  matched_goals = map( g->components_matched(output,g), goallist )
  #println("matched_goals: ",matched_goals)
  fm = findmaxall( map(length, matched_goals ))
  #println("fm: ",fm)
  (Float64(fm[1]), fm[2], matched_goals[fm[2]])
end

# Note that fitness does not take robustness into account
function goals_matched_hamming( output::Vector{MyInt}, goallist::GoalList, numinputs::Int64; 
    avgfitness::Bool=false, perm_heuristic::Bool=false  )
  #println("goals_matched_hamming")
  if perm_heuristic
    goallist_scores = map( g->match_score_perm_heuristic(output,g,numinputs,avgfitness=avgfitness), goallist )
  else
    goallist_scores = map( g->match_score(output,g,numinputs,avgfitness=avgfitness), goallist )
  end
  matched_goals = map( g->components_matched(output,g), goallist )
  #println("goallist_scores: ",goallist_scores)
  #println("matched_goals: ",matched_goals)
  (best_score,best_score_list) = length(goallist_scores)>0 ? findmaxall( goallist_scores ) : (0.0,[])
  (best_score, best_score_list, matched_goals[best_score_list])
end

# Removed robust_sel and active_only keyword args on 6/6/20
# Fault tolerance is the  average Hamming deviation of the output from the unperturbed output under perturbation of the output of each node.
# See equation 3.3 of Macia and Sole (2009)
function next_chromosome!(c::Chromosome, goallist::GoalList, funcs::Vector{Func}, fitness::Float64=0.0; 
      use_robustness::Bool=false, hamming_sel::Bool=true, active_only::Bool=false, num_mutations::Int64=1,
      avgfitness::Bool=false, perm_heuristic=perm_heuristic, fault_tol::Bool=false, ftf_param::Float64=0.95,
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  context = construct_context(c.params.numinputs)
  #println("next_chromosome!:  fault_tol: ",fault_tol,"  ftf_param: ",ftf_param)
  #print_build_chromosome( c )
  if fault_tol 
    old_ftf = fault_tolerance_fitness( c )
    #println("old_ftf: ",old_ftf)
  end
  new_c = deepcopy(c)
  for _ = 1:num_mutations
    mutate_chromosome!( new_c, funcs, insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
  end
  output = execute_chromosome(new_c,context)
  #println("new_c output: ",output,"  new_c.fitness: ",new_c.fitness)
  goals_matched = hamming_sel ? goals_matched_hamming : goals_matched_exact
  ( new_fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs,
      avgfitness=avgfitness, perm_heuristic=perm_heuristic )
  #println("old_fitness: ",fitness,"  c.fitness: ",c.fitness, "  new_fitness: ",new_fitness)
  if !fault_tol
    if new_fitness < fitness
      new_c = c
    else
      new_c.fitness = new_fitness
    end    
  else
    new_ftf = fault_tolerance_fitness( new_c )
    #println("new_ftf: ",new_ftf,"  new_c: ",(new_fitness >= fitness) && (new_ftf >= old_ftf))
    if (new_fitness < fitness) || ((new_ftf < old_ftf) && (rand() > ftf_param))
      new_c = c
    else
      new_c.fitness = new_fitness
    end
  end  
  if use_robustness
    mut_robust = mutational_robustness( new_c, funcs, active_only=active_only )
    new_c.robustness = mut_robust
    #println("new_fitness: ",new_fitness,"  n_mut_robust: ",mut_robust)
  end
  #print_build_chromosome( new_c )
  return ( new_c, matched_goals, matched_goals_list )
end

# Does single-individual neutral evolution.  
# Does subgoal evolution where the number of matched components of a goal is maximized if exact=true
# Does Hamming evolution where the Hamming distance of matched components of a goal is minimized if exact=false
# Removed robust_sel and active_only keyword args on 6/6/20
function mut_evolve( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, max_steps::Integer;
      hamming_sel::Bool=true, use_robustness::Bool=false, num_mutations::Int64=1, print_improvements::Bool=false,
      print_steps::Bool=true,
      avgfitness::Bool=false, perm_heuristic=false, fault_tol::Bool=false, ftf_param::Float64=0.95, 
      fit_limit::Float64=Float64(c.params.numoutputs), 
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  #println("mut_evolve fit limit: ",fit_limit)
  #println("mut evolve avgfitness: ",avgfitness,"  fault_tol: ",fault_tol,"  ftf_param: ",ftf_param)
  #print_build_chromosome(c)
  #print_parameters(c.params)
  #println("length(goallist): ",length(goallist))
  #println("goallist: ",goallist)
  c.fitness = 0.0  # Starting with a non-zero fitness inhibits evolution
  output = output_values(c)   # Executes c if it has not already been executed
  #println("output: ",output)
  #println("trunc(c.fitness): ",trunc(c.fitness))
  goals_matched = hamming_sel ? goals_matched_hamming : goals_matched_exact
  ( fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs, avgfitness=avgfitness )
  fitness = c.fitness
  orig_c = deepcopy(c)
  #println("initial fitness: ",fitness)
  robustness = 0.0
  step = 0
  worse = 0
  same = 0
  better = 0
  prev_c = deepcopy(c)  # previous generation c
  @assert output_values(prev_c) == output_values(orig_c)
  #(c, new_robustness, matched_goals_list ) = 
  #    next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, avgfitness=avgfitness ) 
  (c, matched_goals, matched_goals_list ) = 
      next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, 
      perm_heuristic=perm_heuristic, avgfitness=avgfitness, fault_tol=fault_tol, ftf_param=ftf_param, insert_gate_prob=insert_gate_prob, 
      delete_gate_prob=delete_gate_prob ) 
  output = output_values(c)   # Executes c if it has not already been executed
  #println("output: ",output)
  #while step < max_steps && new_fitness < c.params.numoutputs
  #println("c.fitness: ",c.fitness,"  fit_limit: ",fit_limit)
  while step < max_steps && c.fitness < fit_limit
    #println("step: ",step,"  output: ",output,"   ")
    #println("trunc(c.fitness): ",trunc(c.fitness))
    if c.fitness > fitness
      if print_improvements
        matched_goal_components = [map( x->x[1], mgl) for mgl in matched_goals_list]
        println("fitness improved from ",fitness," to ",c.fitness," at step: ",step, "  goals_matched: ", matched_goals, "  matched goal componets: ", matched_goal_components)
      end
      fitness = c.fitness 
      worse = 0
      same = 0
      better += 1
    elseif c.fitness == fitness 
      #println("new chromosome with fitness: ",fitness )
      same += 1
    else
      #println("discarded chromosome with new_fitness: ",c.fitness)
      c = prev_c
      worse += 1
    end
    if prev_c.fitness < orig_c.fitness
      #println("step: ",step,"  prev_c.fitness: ",prev_c.fitness,"  orig_c.fitness: ",orig_c.fitness)
    end
    step += 1
    prev_c = deepcopy(c)
    if step < max_steps   # Construct c for next generation
      (c, matched_goals, matched_goals_list ) = 
          next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, 
          perm_heuristic=perm_heuristic, avgfitness=avgfitness, fault_tol=fault_tol, ftf_param=ftf_param, 
          insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob ) 
            #next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, use_robustness=use_robustness, 
            #num_mutations=num_mutations, perm_heuristic=perm_heuristic, avgfitness=avgfitness, fault_tol=fault_tol,
            #ftf_param=ftf_param ) 
      output = output_values(c)   # Executes c if it has not already been executed
    end
  end
  if step == max_steps   # Failed to find goal
    if print_steps
      println("mut_evolve finished at step limit ",max_steps," with fitness: ", c.fitness ) 
    end
  else
    matched_g = map( x->x[1][1], matched_goals_list )
    if print_steps
      println("mut_evolve finished in ",step," steps for goal ",matched_g," with fitness: ", c.fitness )
    end
    #println("matched_goals: ",matched_goals,"  matched_goals_list: ",matched_goals_list)
  end
  if orig_c.fitness > c.fitness   # this should never happen
    #println("(orig_c.fitness,c.fitness): ",(orig_c.fitness,c.fitness)) 
    output = output_values(orig_c)
    ( fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs )
    return(orig_c,step,0,0,0,output,matched_goals,matched_goals_list)
  end
  ##println("worse: ",worse,"  same: ",same,"  better: ",better)
  #println("matched_goals: ",matched_goals)
  if step < max_steps
    if print_improvements
      matched_goal_components = [map( x->x[1], mgl) for mgl in matched_goals_list]
      println("fitness improved from ",fitness," to ",c.fitness," at step: ",step, "  goals_matched: ", matched_goals, "  matched goal comps: ", matched_goal_components)
    end
    #println("sort(output: ",sort(output),"  sort(goallist[matched_goals[1]]): ",sort(goallist[matched_goals[1]]))
    if Int(trunc(fit_limit)) == c.params.numoutputs
      @assert sort(output) == sort(goallist[matched_goals[1]])
    end
  end
  (c,step,worse,same,better,output,matched_goals,matched_goals_list)
end 

# Does up to n_repeats tries to evolve a chromosome which maps to a goal in goallist
function mut_evolve_repeat(n_repeats::Int64, p::Parameters, goallist::GoalList, funcs::Vector{Func}, max_steps::Integer;
      hamming_sel::Bool=true, use_robustness::Bool=false, num_mutations::Int64=1, print_improvements::Bool=false,
      avgfitness::Bool=false, perm_heuristic=false, fault_tol::Bool=false, ftf_param::Float64=0.95, 
      fit_limit::Float64=Float64(p.numoutputs) )
  repeat = 1
  c = random_chromosome( p, funcs )
  (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, goallist, funcs, max_steps, 
        hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, print_improvements=print_improvements,
        avgfitness=avgfitness, perm_heuristic=perm_heuristic, fault_tol=fault_tol, ftf_param=ftf_param,
        fit_limit=fit_limit ) 
  while repeat < n_repeats && step == max_steps
    c = random_chromosome( p, funcs )
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, goallist, funcs, max_steps, 
        hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, print_improvements=print_improvements,
        avgfitness=avgfitness, perm_heuristic=perm_heuristic, fault_tol=fault_tol, ftf_param=ftf_param,
        fit_limit=fit_limit ) 
    repeat += 1
  end
  if step < max_steps  # circut that outputs a goal of goallist has been found
    return (c,step,worse,same,better,output,matched_goals,matched_goals_list)
  else
    return nothing
  end
end

# If mut_evolve() fails in max_steps iterations, increase the numinteriors and numlevelsback and try again.
# Continue up to repeat_limit until it suceeds.
function mut_evolve_increase_numints( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, max_steps::Integer;
      hamming_sel::Bool=true, use_robustness::Bool=false, num_mutations::Int64=1, print_improvements::Bool=false,
      avgfitness::Bool=false, perm_heuristic=false, fault_tol::Bool=false, ftf_param::Float64=0.95,
      fit_limit::Float64=Float64(c.params.numoutputs) ) 
  total_steps = 0
  new_numints = c.params.numinteriors
  new_levsback = c.params.numlevelsback
  repeat_limit = 10  
  numinteriors_repeat = c.params.numinteriors
  numlevelsback_repeat = c.params.numlevelsback
  (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, goallist, funcs, max_steps, 
      hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, print_improvements=print_improvements,
      avgfitness=avgfitness, perm_heuristic=perm_heuristic, fault_tol=fault_tol, ftf_param=ftf_param,
      fit_limit=fit_limit ) 
  repeat = 0
  while repeat < repeat_limit && step == max_steps
    #println("starting repeat loop repeat: ",repeat)
    total_steps += max_steps
    numinteriors_repeat += 1
    numlevelsback_repeat += 1
    p_repeat = Parameters( numinputs=c.params.numinputs, numoutputs=c.params.numoutputs, numinteriors=numinteriors_repeat,
        numlevelsback=numlevelsback_repeat )
    println("repeating function mut_evolve with numints: ",numinteriors_repeat,"  and with levsback: ",numlevelsback_repeat)
    c = random_chromosome( p_repeat, funcs )
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, goallist, funcs, max_steps, 
        hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, print_improvements=print_improvements,
        avgfitness=avgfitness, perm_heuristic=perm_heuristic, fault_tol=fault_tol, ftf_param=ftf_param,
        fit_limit=fit_limit ) 
    new_numints = p_repeat.numinteriors
    new_levsback = p_repeat.numlevelsback
    repeat += 1
  end
  total_steps += step
  if repeat == repeat_limit
    error(" repeat reached repeat_limit in function mut_evolve_increase_numints()")
  end
  (c,total_steps,worse,same,better,output,matched_goals,matched_goals_list,new_numints,new_levsback)
end

# Attempts to reduce the number active of the chromosome c by mutational neutral evolution
# The chomosome c is assumed to be at optimal fitness, and this function will terminate with
#   an error if it is not optimal.
# Then mutation is done nreps times, and the mutation is accepted if number_active(c) is reduced.
function mut_reduce_numactive( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, nreps::Integer;
      hamming_sel::Bool=true, num_mutations::Int64=1, print_improvements::Bool=false )
  output = output_values(c)   # Executes c if it has not already been executed
  goals_matched = hamming_sel ? goals_matched_hamming : goals_matched_exact
  ( fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs )
  if sort(output) != sort(goallist[matched_goals[1]])
    error("input chromosome must be have optimal fitness in function mut_reduce_numactive()")
  end
  best_num_active = number_active(c)
  #println("init num_active  num_active: ",best_num_active)
  orig_c = deepcopy(c)
  #println("initial fitness: ",fitness)
  step = 0
  worse = 0
  same = 0
  better = 0
  prev_c = deepcopy(c)  # previous generation c
  best_c = c
  @assert output_values(prev_c) == output_values(orig_c)
  (c, matched_goals, matched_goals_list ) = next_chromosome!(c, goallist, funcs, 
      hamming_sel=hamming_sel, use_robustness=false, num_mutations=num_mutations ) 
  new_num_active = number_active(c)
  #println("before while new_num_active: ",new_num_active)
  #while step < nreps && new_fitness < c.params.numoutputs
  while step < nreps 
    #println("step: ",step,"  output: ",output,"   ")
    if c.fitness > fitness  # This should never happen
      error("c.fitness > fitness")
      if print_improvements
        matched_goal_components = [map( x->x[1], mgl) for mgl in matched_goals_list]
        println("fitness improved from ",fitness," to ",c.fitness," at step: ",step, "  goals_matched: ", matched_goals, "  matched goal components: ", matched_goal_components)
      end
      fitness = c.fitness 
      worse = 0
      same = 0
      better += 1
    elseif c.fitness == fitness && new_num_active < best_num_active
      #println("new chromosome with fitness: ",fitness, "  and num_active: ", new_num_active )
      if print_improvements
        matched_goal_components = [map( x->x[1], mgl) for mgl in matched_goals_list]
        println("num_active decreased from ",best_num_active," to ",new_num_active," at step: ",step,"  goals_matched: ", matched_goals, "  matched goal componets: ", matched_goal_components)
      end
      worse = 0
      same = 0
      better += 1
      best_c = deepcopy(c)
      @assert new_num_active == number_active(c)
      best_num_active = new_num_active
      #println("num_active_improved  new_num_active: ",new_num_active)
    elseif c.fitness == fitness && new_num_active == best_num_active
      #println("new chromosome with fitness: ",fitness, "  and num_active: ", new_num_active )
      same += 1
    elseif c.fitness == fitness && new_num_active > best_num_active
      #println("discarded chromosome with fitness: ",fitness, " and num_active: ", new_num_active )
      c = prev_c
      worse += 1
    elseif c.fitness < fitness
      #println("discarded chromosome with new_fitness: ",c.fitness)
      c = prev_c
      worse += 1
    end
    step += 1
    prev_c = deepcopy(c)
    if step < nreps
      (c, matched_goals, matched_goals_list ) = next_chromosome!(c, goallist, funcs, 
            hamming_sel=hamming_sel, use_robustness=false, num_mutations=num_mutations ) 
      #output = output_values(c)   # Executes c if it has not already been executed
      new_num_active = number_active(c)
      #println("end while new_num_active: ",new_num_active)
    end
  end
  if step == nreps
    #println("mut_reduce_numactive finished at step limit ",nreps," with fitness: ", best_c.fitness, " and num_active: ", best_num_active ) 
  else
    error("step == nreps")
    println("mut_reduce_numactive finished in ",step," steps with fitness: ", c.fitness, " and num_active: ", new_num_active )
  end
  if orig_c.fitness > c.fitness
    error("orig_c.fitness > c.fitness")
    #println("(orig_c.fitness,c.fitness): ",(orig_c.fitness,c.fitness)) 
    output = output_values(orig_c)
    ( fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs )
    return(orig_c,step,0,0,0,output,matched_goals,matched_goals_list)
  end
  ##println("worse: ",worse,"  same: ",same,"  better: ",better)
  #println("matched_goals: ",matched_goals)
  if step < nreps
    if print_improvements
      matched_goal_components = [map( x->x[1], mgl) for mgl in matched_goals_list]
      println("num_active decreased from ",num_active," to ",new_num_active," at step: ",step, "  goals_matched: ", matched_goals, "  matched goal componets: ", matched_goal_components)
    end
    @assert sort(output) == sort(goallist[matched_goals[1]])
  end
  ( fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, best_c.params.numinputs )
  #println("assert  best_num_active: ",best_num_active,"  number_active(best_c): ",number_active(best_c))
  @assert best_num_active == number_active(best_c)
  (best_c,step,worse,same,better,output,matched_goals,matched_goals_list)
end 

# Do a random walk through circuit space starting with chromosome c which must be optimal for the given goallist.
# This function is trying to minimize num_active.  
# There is another random_neutral_walk function in Evolability.jl which is trying to measure evolability
function random_neutral_walk( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, nreps::Integer;
      hamming_sel::Bool=true, num_mutations::Int64=1, print_improvements::Bool=false, reduce_numactive_reps::Int64=0 )
  sum_num_active = 0
  min_num_active = 2^20  # a very large integer
  count_saved_reps = 0
  count_min_num_active = 0
  min_active_chromes_list = Chromosome[]
  output = output_values(c)   # Executes c if it has not already been executed
  goals_matched = hamming_sel ? goals_matched_hamming : goals_matched_exact
  ( orig_fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs )
  if sort(output) != sort(goallist[matched_goals[1]])
    error("input chromosome must be have optimal fitness in function mut_reduce_numactive()")
  end
  println("orig_fitness: ",orig_fitness)
  prev_c = deepcopy(c)
  for i = 1:nreps
    (c, matched_goals, matched_goals_list ) = next_chromosome!(c, goallist, funcs, 
          hamming_sel=hamming_sel, use_robustness=false, num_mutations=num_mutations ) 
    if c.fitness < orig_fitness 
      #println("discarding chromosome at rep: ",i)
      c = prev_c
      continue
    end
    new_num_active = number_active(c)
    count_saved_reps += 1
    #println("rep: ",i,"  num_active: ",new_num_active)
    sum_num_active += new_num_active
    if min_num_active > new_num_active
      println("new min: ",new_num_active,"  at rep: ",i)
      min_num_active = new_num_active
      count_min_num_active = 1
      @assert sort(output_values(c)) == sort(goallist[matched_goals[1]])
      min_active_chromes_list = Chromosome[c]
      if reduce_numactive_reps > 0
        (best_c,step,worse,same,better,output,matched_goals,matched_goals_list) =
            mut_reduce_numactive( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, nreps::Integer, num_mutations=num_mutations )
        best_num_active = number_active(best_c)
        if best_num_active < new_num_active
          println("best_num_active: ",best_num_active)
        end
      end
    elseif min_num_active == new_num_active
      #println("alt min: ",new_num_active,"  at rep: ",i)
      @assert sort(output_values(c)) == sort(goallist[matched_goals[1]])
      count_min_num_active += 1
      if length(min_active_chromes_list) < 4
        Base.push!(  min_active_chromes_list, c )
      end
      if reduce_numactive_reps > 0 && count_min_num_active < 10
        (best_c,step,worse,same,better,output,matched_goals,matched_goals_list) =
            mut_reduce_numactive( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, nreps::Integer, num_mutations=num_mutations )
        best_num_active = number_active(best_c)
        if best_num_active < new_num_active
          println("best_num_active: ",best_num_active)
        end
      end
    end
    prev_c = deepcopy(c)
  end
  println("nreps: ",nreps,"  nmuts: ",num_mutations,"  avg num active: ",sum_num_active/count_saved_reps,
        "  min_num_active: ",min_num_active, "  count_min_num_active: ",count_min_num_active)
  return (min_num_active,count_min_num_active,min_active_chromes_list)
end

# Given a goal, for each iteraion of numiterations,
#    evolve two chromosomes that map to it, and try to find a distance nondecreasing path from one to the other
# TODO:  What happens when mut_evolve to find c1 or c2 fails?  Retry?  How many times?
function run_circuit_evolve( goal::Goal, p::Parameters, numiterations::Int64, maxreps::Int64 )
  funcs=default_funcs(p.numinputs)
  successes = 0
  sum_steps = 0
  for iter = 1:numiterations
    c = random_chromosome(p,funcs)
    res = mut_evolve(c,[goal],funcs,maxreps); c1=res[1]
    c = random_chromosome(p,funcs)
    res = mut_evolve(c,[goal],funcs,maxreps); c2=res[1]
    steps = circuit_evolve(c1,c2, maxreps )
    if steps <= maxreps
      println("success with ",steps," steps")
      sum_steps += steps
      successes += 1
    end
  end
  (successes,sum_steps/numiterations)
end
    
# Do neutral evolution starting at chromosome c_src trying to get to chromosome c_dest.
# c_src and c_dest must have the same parameters and the same output values (goals).
function circuit_evolve( c_src::Chromosome, c_dest::Chromosome, maxreps::Int64 )
  funcs = default_funcs(c_src.params.numinputs)
  @assert c_src.params == c_dest.params
  outputs_src = output_values( c_src )
  outputs_dest = output_values( c_dest )
  @assert outputs_dest == outputs_src
  c = deepcopy(c_src)  # current chromosome
  dist = circuit_distance( c_src, c_dest )
  #outputs_prev = output_values( c_src )
  step = 1
  c_new = deepcopy(c)
  mutate_chromosome!( c_new, funcs )
  #print("c_new before loop: ")
  #print_build_chromosome(c_new)
  while step <= maxreps
    if outputs_src != output_values(c_new)
      step += 1
      #println("continue")
      #continue  # continue with current ciruit c unchanged
    else
      println("c_new code:  ",circuit_code(c_new))
      println("c_dest_code: ",circuit_code(c_dest))
      new_dist = circuit_distance( c_new, c_dest )
      println("dist:     ",dist,"  new_dist: ",new_dist)
      if new_dist > dist
        step += 1
        #continue  # continue with current ciruit c unchanged
      elseif new_dist == 0  # success
        println("circuit_evolve() found a path to c_dest!")
        return step
      elseif new_dist < dist
        println("distance improved from ",dist," to ",new_dist)
        dist = new_dist
        c = deepcopy(c_new)
        step += 1
      end
    end
    c_new = deepcopy(c)
    mutate_chromosome!(c_new,funcs)
    #print("c_new end of loop: ")
    #print_build_chromosome(c_new)
  end
  println("circuit_evolve() failed to find a path to c_dest!")
  return step
end
  
# Change the fields of chromosome c to be the fields of chromosom c_to_copy
# TODO:  Move to Chromosome.jl
function copy_chromosome!( c::Chromosome, c_to_copy::Chromosome )
  c.params = c_to_copy.params
  c.inputs = c_to_copy.inputs
  c.interiors = c_to_copy.interiors
  c.outputs = c_to_copy.outputs
  c.interiors = c_to_copy.interiors
end

#=  Commented out on 4/17/21
# Example call:  build_chromosome((1,2), ((OR,[1,2]),(AND,[2,3])),(4,),0.0)
function build_chromosome( inputs::Tuple, ints::Tuple, outs::Tuple, fitness::Float64 )
  num_in = length(inputs)
  num_ints = length(ints)
  num_outs = length(outs)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=num_ints+num_outs )
  in_nodes = [InputNode(in_index) for in_index in inputs]
  int_nodes = [InteriorNode(int_pair[1], int_pair[2]) for int_pair in ints]
  out_nodes = [OutputNode(out_index) for out_index in outs]
  Chromosome( p, in_nodes, int_nodes, out_nodes, fitness, 0.0 )
end
=#

# The next functions should be in a Utilities.jl file
#  Returns a double whose first component is the maximum value in X, and whose second component is
#     the list of indices that give the maximum.
#=
function findmaxall( X::AbstractVector )
  max = findmax(X)[1]
  #println("max: ",max)
  Xind = [(X[i],i) for i = 1:length(X)]
  #println("Xind:",Xind)
  Xmax = filter(x->(x[1]==max),Xind)
  #println("Xmax:",Xmax)
  indices = [x[2] for x in Xmax]
  (max,indices)
end
=#

# argument A is a vector where > and == can be used to compare elements
function findmaxall( A::AbstractVector )
  if isempty(A)
    error("Empty list in findmaxall()")
  end
  indices = [1]
  for i = 2:length(A) 
    if A[i] > A[indices[1]]
      indices = [i]
    elseif A[i] == A[indices[1]]
      push!(indices,i)
    end
  end
  (A[indices[1]],indices)
end
 
# argument A is a vector where < and == can be used to compare elements
function findminall( A::AbstractVector )
  if isempty(A)
    error("Empty list in findminall()")
  end
  indices = [1]
  for i = 2:length(A) 
    if A[i] < A[indices[1]]
      indices = [i]
    elseif A[i] == A[indices[1]]
      push!(indices,i)
    end
  end
  (A[indices[1]],indices)
end

# cmp is a two-argument comparison function:  greater-than for max, less-than for min
# eq  is a two-argument equality function
function findmaxminall( cmp, eq, A::AbstractVector )
  if isempty(A)
    error("Empty list in findmaxminall()")
  end
  indices = [1]
  for i = 2:length(A) 
    if cmp(A[i], A[indices[1]])
      indices = [i]
    elseif eq(A[i], A[indices[1]])
      push!(indices,i)
    end
  end
  (A[indices[1]],indices)
end

# findmaxall where the comparison functions look at the first components of the elements of A
findmaxall_1( A::AbstractVector ) = findmaxminall( (x,y)->(x[1]>y[1]), (x,y)->(x[1]==y[1]), A )

# findmnall where the comparison functions look at the first components of the elements of A
findminall_1( A::AbstractVector ) = findmaxminall( (x,y)->(x[1]<y[1]), (x,y)->(x[1]==y[1]), A )
 

#  Returns a double whose first component is the minimum value in X, and whose second component is
#     the list of indices that give the minimum.
function findminall( X::AbstractVector )
  min = findmin(X)[1]
  #println("min: ",min)
  Xind = [(X[i],i) for i = 1:length(X)]
  #println("Xind:",Xind)
  Xmin = filter(x->(x[1]==min),Xind)
  #println("Xmin:",Xmin)
  indices = [x[2] for x in Xmin]
  (min,indices)
end

# Returns the maximum value in A along with a random index for that value
function findmaxrand( A::AbstractVector )  
  (val,indices) = findmaxall( A )
  if length(indices) == 1
    return (val,indices[1])
  else
    ind = rand(indices)
    #println("indices: ",indices,"  ind: ",ind)
    return (val,ind)
  end
end

# Returns the minimum value in A along with a random index for that value
function findminrand( A::AbstractVector )  
  (val,indices) = findmaxall( A )
  if length(indices) == 1
    return (val,indices[1])
  else
    ind = rand(indices)
    #println("indices: ",indices,"  ind: ",ind)
    return (val,ind)
  end
end

# Use neutral evolution to find a chromosome c that maximizes funct(c)
# funct( c::Chromosome ) returns a Float64.
function evolve_function( funct::Function, p::Parameters, funcs::Vector{Func}, max_steps::Int64;
    goallist::Vector{Vector{MyInt}}=[MyInt[]], max_evolve_steps::Int64=10000 )
  println("goallist: ",goallist)
  use_goal = !(goallist==[MyInt[]])
  c = random_chromosome( p, funcs )
  orig_c = deepcopy(c)
  goal = goallist[1]
  println("goal: ",goal)
  if use_goal
    result = mut_evolve( c, goallist, funcs, max_evolve_steps, print_improvements = true ) 
    steps = result[2]
    if steps == max_evolve_steps
      error("evolution to goal failed in function evolve_function()")
    else
      println("evolution to goal suceeded in function evolve_function()")
    end 
    c = result[1]
    goal = output_values(c)
    orig_c = deepcopy(c)
  end
  current_fitness = funct(c)
  println("starting fitness: ",current_fitness)
  for i = 1:max_steps
    prev_c = deepcopy(c)
    (c,active) = mutate_chromosome!(c,funcs)
    if use_goal && output_values(c) != goal
      c = prev_c
      continue
    end  
    new_fitness = funct(c)
    if new_fitness < current_fitness 
      c = prev_c
    elseif new_fitness > current_fitness
      println("i: ",i,"  new_fitness: ", new_fitness,"  funct(c): ",funct(c) )
      current_fitness = new_fitness
    end
    println("i: ",i,"  current_fitness: ", current_fitness, "  funct(c): ",funct(c)  )
  end
  (c,orig_c)
end

# to test circuit_evolve()
function test_evolve()
  p = Parameters(numinputs=2,numoutputs=1,numinteriors=3,numlevelsback=3); 
  funcs=default_funcs(p.numinputs);c = random_chromosome(p,funcs); goal=output_values(c)
  res = mut_evolve(c,[goal],funcs,1000); c_dest=res[1]
  circuit_evolve(c,c_dest, 1000 )
end

# Evolves a chromosome (cirucit) that maps to g starting with chromosome c.
# max_steps is the maximum number of evolutionary steps.
# If evolution hasn't succeeeded in max_steps, return nothing.
# insert_gate_prob is the probability of inserting a gate on a mutation of a chromosome.
# delete_gate_prob is similar for deleting a gate.
# Similar to mut_evolve except that this takes a single goal instead of a goal list as an argument.
function neutral_evolution( c::Chromosome, g::Goal, max_steps::Integer; print_steps::Bool=false,
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  funcs = default_funcs( c.params.numinputs )
  step = 0
  ov = output_values( c) 
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  new_c = deepcopy(c)
  while step < max_steps && ov != g
    step += 1
    (new_c,active) = mutate_chromosome!( new_c, funcs, insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
    new_ov = output_values( new_c )
    new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    #println("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur dis: ",current_distance,"  new_dis: ",new_distance )
    if new_ov == ov 
      c = new_c
      if print_steps
        println("step: ",step," is neutral.")
      end
    elseif new_distance < current_distance
      if print_steps
        println("step: ",step,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance)
      end
      c = new_c
      ov = new_ov
      current_distance = new_distance
    else
      if print_steps
        print("step: ",step,"  new_output: ",new_ov,"  new circuit: ")
        print_circuit( new_c )
      end 
    end
  end
  if step == max_steps
    #println("neutral evolution failed with ",step," steps for goal: ",g)
    return (nothing, step)
  else
    #println("neutral evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end

# Evolves a LinCirucit that maps to g starting with chromosome c.
# max_steps is the maximum number of evolutionary steps.
# If evolution hasn't succeeeded in max_steps, return nothing.
# insert_gate_prob is the probability of inserting a gate on a mutation of a chromosome.
# delete_gate_prob is similar for deleting a gate.
# Similar to mut_evolve except that this takes a single goal instead of a goal list as an argument.
function neutral_evolution( c::Circuit, g::Goal, max_steps::Integer; print_steps::Bool=false,
      insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0 )
  LinCirc = typeof(c) == LinCircuit ? :true : :false
  funcs = default_funcs( c.params.numinputs )
  step = 0
  ov = output_values( c) 
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  new_c = deepcopy(c)
  while step < max_steps && ov != g
    step += 1
    if typeof(c) == Chromosome
      (new_c,active) = mutate_chromosome!( new_c, funcs, insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
    elseif typeof(c) == LinCircuit
      new_c = mutate_circuit!( new_c, funcs )
    end
    new_ov = output_values( new_c )
    new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    #println("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur dis: ",current_distance,"  new_dis: ",new_distance )
    if new_ov == ov 
      c = new_c
      if print_steps
        println("step: ",step," is neutral.")
      end
    elseif new_distance < current_distance
      if print_steps
        println("step: ",step,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance)
      end
      c = new_c
      ov = new_ov
      current_distance = new_distance
    else
      if print_steps
        print("step: ",step,"  new_output: ",new_ov,"  new circuit: ")
        if LinCirc
          print_circuit( new_c, funcs )
        else
          print_circuit( new_c )
        end
      end 
    end
  end
  if step == max_steps
    #println("neutral evolution failed with ",step," steps for goal: ",g)
    return (nothing, step)
  else
    #println("neutral evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end

# Computes phenotype evolvability, genotype evolvability, robustness, complexity, steps
# num_circuits is the number of circuits used to compute properties.
function geno_circuits( g::Goal, p::Parameters, num_circuits::Integer, max_steps::Integer, max_attempts::Integer )
  funcs = default_funcs( p.numinputs )
  c = random_chromosome( p, funcs )
  sum_steps = 0
  circuit_list = Chromosome[]
  n_circuits = 0
  attempt = 0
  nc = nothing
  while attempt < max_attempts && n_circuits < num_circuits
    attempt += 1
    (nc,steps) = neutral_evolution( c, g, max_steps )
    sum_steps += steps
    if nc != nothing
      n_circuits += 1
      push!( circuit_list, nc )
    end
  end
  if n_circuits < num_circuits
    println("geno_properties failed to find num_circuits circuits mapping to goal: ",g," in ", attempt," attempts.")
    return (nothing, sum_steps)
  end
  return (circuit_list, sum_steps)
end
 
function geno_properties( cl_ss::Tuple{Vector{Chromosome}, Int64} ) # Tuple of circuit list and sum_steps.
  (cl, sum_steps) = cl_ss
  funcs = default_funcs( cl[1].params.numinputs )
  g = output_values(cl[1])
  sum_robust = 0.0
  sum_complexity = 0.0
  sum_evolvability = 0.0  # genotype evolability
  genotypes = Goal[]
  genotype_set = Set(Goal[])
  n_circuits = length(cl)
  for c in cl
    @assert output_values(c) == g
    sum_robust += mutational_robustness( c, funcs )
    sum_complexity += complexity5(c) 
    genotypes = Set(mutate_all(c,funcs,output_outputs=true,output_chromosomes=false))
    sum_evolvability += length(genotypes)-1
    #println("len(genotypes): ",length(genotypes),"  genotypes: ",genotypes)
    genotype_set = union( genotype_set, genotypes )
  end
  (g, sum_robust/n_circuits, sum_complexity/n_circuits, sum_evolvability/n_circuits, length(genotype_set), sum_steps/n_circuits) 
end

# Computes the properites robustness, complexity, genotype evolability, phenotype evolvability, evolutonary steps
#    for the goals of gl.  The properties are based on the evoution of num_circuits circuits that map to the goal.
# max_attempts is the maximum number of calls to neutral_evolution() that are made in the attempt to evolve num_circuits circuits
#    that evolve to a goal 
# max_steps is the maximum number of evolutionary steps done by neutral_evoluion() during 1 attempt to evolve a circuit
#    that maps to a goal.
function geno_list_properties( gl::GoalList, p::Parameters, num_circuits::Integer, max_steps::Integer, max_attempts::Integer;
    csvfile::String="" )
  df = DataFrame()
  df.goal = Goal[]
  df.robustness = Float64[]
  df.complexity = Float64[]
  df.g_evolvability = Float64[]
  df.p_evolvability = Float64[]
  df.steps = Float64[]
  gp_list = Vector{Tuple{Goal,Float64,Float64,Float64,Float64,Float64}}[]
  gp_list = pmap( g->geno_properties( geno_circuits( g, p, num_circuits, max_steps, max_attempts ) ), gl )
  for gp in gp_list
    push!(df, gp) 
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# num_circuits: ",num_circuits)
      println(f,"# max_steps: ",max_steps)
      println(f,"# max_attempts: ",max_attempts)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
