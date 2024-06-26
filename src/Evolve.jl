# Restored to stash/Evolve12_15_22.jl on 12/19/22.  Gave up on adding robustness:  too many failures.
#using .Threads
export components_matched, goals_matched_exact, goals_matched_hamming, next_chromosome!, mut_evolve, mut_evolve_increase_numints
export mut_reduce_numactive
export random_neutral_walk, match_score, findmaxall, findminall, findmaxall_1, findminall_1, findmaxrand, findminrand
export evolve_function, mut_evolve_repeat, circuit_evolve, run_circuit_evolve, run_ph_evolve
export neutral_evolution, geno_circuits, geno_properties, geno_list_properties, lambda_evolution, pheno_evolve, run_pheno_evolve
export run_evolve_to_pheno_mt, run_to_rand_phenos_mt
export directed_neutral_evolution
export epochal_evolution_fitness, run_epochal_evolution_fitness
export rand_bit_word

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
    #print("next chromosome: ")
    #print_circuit(new_c)
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
  default_funcs(c.params)   # Set the global variable Ones when running on multiple processes.
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
  #print_circuit(c)
  output = output_values(c)   # Executes c if it has not already been executed
  #println("output: ",output)
  #while step < max_steps && new_fitness < c.params.numoutputs
  #println("c.fitness: ",c.fitness,"  fit_limit: ",fit_limit)
  while step < max_steps && c.fitness < fit_limit
    #print("step: ",step,"  output: ",output," ")
    #print_circuit(c)
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
      #print_circuit(c)
      worse += 1
    end
    if prev_c.fitness < orig_c.fitness
      println("step: ",step,"  prev_c.fitness: ",prev_c.fitness,"  orig_c.fitness: ",orig_c.fitness)
    end
    step += 1
    prev_c = deepcopy(c)
    #println("end loop ",rand(1:100))
    if step < max_steps   # Construct c for next generation
      (c, matched_goals, matched_goals_list ) = 
          next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, 
            perm_heuristic=perm_heuristic, avgfitness=avgfitness, fault_tol=fault_tol, ftf_param=ftf_param, 
            insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob ) 
            #next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, use_robustness=use_robustness, 
            #num_mutations=num_mutations, perm_heuristic=perm_heuristic, avgfitness=avgfitness, fault_tol=fault_tol,
            #ftf_param=ftf_param ) 
      output = output_values(c)   # Executes c if it has not already been executed
      #println("step: ",step+1,"  output: ",output)
      #print_circuit(c)
    end
  end
  if step == max_steps   # Failed to find goal
    if print_steps
      println("mut_evolve finished at step limit ",max_steps," with fitness: ", c.fitness ) 
    end
  else
    matched_g = map( x->[x[i][1] for i = 1:c.params.numoutputs], matched_goals_list )
    if print_steps
      println("mut_evolve finished in ",step," steps for goal ",matched_g[1]," with fitness: ", c.fitness )
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
function mut_evolve_increase_numints( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, max_steps::Integer, n_repeats::Int64=10;
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
    println("repeating function mut_evolve_increse numints: ",numinteriors_repeat,"  and with levsback: ",numlevelsback_repeat,
        "  goallist[1]: ",goallist[1])
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

# Run pheno_evolve() to evolve numcircuits circuits that map to each phenotype in phlist.
# Does max_tries attempts to evolve one circuit that maps to each phenotype in phlist
# Not currently tested for LGP.
function run_ph_evolve( p::Parameters, funcs::Vector{Func}, phlist::GoalList, numcircuits::Int64, max_tries::Int64, max_steps::Int64; 
    use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false, csvfile::String="" )
  if max_tries <= numcircuits
    error("function run_ph_evolve(): max_tries should be greater than numcircuits. Your values: max_tries: ",max_tries,"  numcircuits: ",numcircuits)
  end
  println("methods(robustness): ",methods(robustness))
  rdict = redundancy_dict(p,funcs)
  kdict = kolmogorov_complexity_dict(p,funcs)
  k_comp_funct( x::MyInt ) = typeof(kdict) <: Dict ? kdict[x] : 0
  redund_funct( x::MyInt ) = typeof(rdict) <: Dict ? rdict[x] : 0
  nphenos = length(phlist)
  #result_list = map( ph->pheno_evolve( p, funcs, ph, numcircuits, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), phlist ) 
  result_list = pmap( ph->pheno_evolve( p, funcs, ph, numcircuits, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), phlist ) 
  mean_steps_list = Float64[]
  median_steps_list = Float64[]
  max_steps_list = Float64[]
  std_steps_list = Float64[]
  robustness_list = Float64[]
  log_redundancy_list = Float64[]
  K_complexity_list = Int64[]
  #println("result_list: ",result_list)
  for i = 1:length(result_list)
    steps_list = map(x->x[2], result_list[i][1:numcircuits])
    #println("steps_list: ",steps_list)
    mean_steps = mean(steps_list)
    push!(mean_steps_list,mean_steps)
    median_steps = median(steps_list)
    push!(median_steps_list,median_steps)
    max_steps = findmax(steps_list)[1]
    push!(max_steps_list,mean_steps)
    std_steps = std(steps_list)
    push!(robustness_list,mean( robustness( result_list[i][j][1], funcs ) for j = 1:numcircuits ) )
    push!(std_steps_list,std_steps)
    push!(log_redundancy_list,lg10(redund_funct(phlist[i][1]) ))
    push!(K_complexity_list,k_comp_funct(phlist[i][1]) )
  end
  mean_steps_list,
  median_steps_list,
  std_steps_list
  df = DataFrame( :phlist=>phlist, :numinputs=>fill(p.numinputs,nphenos), :numgates=>fill(p.numinteriors,nphenos), :levsback=>fill(p.numlevelsback,nphenos), 
      :mean_steps=>mean_steps_list, :median_steps=>median_steps_list, :max_steps=>max_steps_list, :std_steps=>std_steps_list, :Kcomp=>K_complexity_list, 
      :robustness=>robustness_list, :lg_redund=>log_redundancy_list )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# MyInt: ",MyInt)
      println(f,"# ",use_lincircuit ? "LGP" : "CGP" )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# length(phlist): ",length(phlist))
      println(f,"# numcircuits: ",numcircuits)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

#run_ph_evolve( p, funcs, map(x->[x],collect(MyInt(0):MyInt(2^2^p.numinputs-1))), 3, 5, 100_000 )

# Run pheno_evolve() to evolve one circuit that maps to each phenotype in phlist.
# Does max_tries attempts to evolve one circuit that maps to each phenotype in phlist
# Not currently tested for LGP.
function run_pheno_evolve( p::Parameters, funcs::Vector{Func}, phlist::GoalList, max_tries::Int64, max_steps::Int64; 
    use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false, csvfile::String="" ) 
  function ph_evolve( p::Parameters, funcs::Vector{Func}, goal::Goal, max_tries::Int64, max_steps::Int64; 
      use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false ) 
    # This is the version of pheno_evolve that finds one circuit per goal.
    (c,steps,total_failures) = pheno_evolve( p, funcs, goal, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps )
    return c != nothing ? (c,steps,number_active(c),robustness(c,funcs),total_failures) : (nothing,nothing,nothing,nothing,nothing)
  end
  println("length(funcs): ",length(funcs))
  df = DataFrame( :numinputs=>Int64[], :numgates=>Int64[], :numlevelsback=>Int64[], :fail_fract=>Float64[], :mean_steps=>Float64[], :median_steps=>Float64[], 
      :std_steps=>Float64[], :mean_nactive=>Float64[], :mean_robust=>Float64[], :first_fail=>Int64[], :subseqent_fail=>Int64[], :total_failures_list=>Vector{Vector{Int64}}[] )
  #result_list = filter(x->x[1]!=nothing, Folds.map( ph->ph_evolve( p, funcs, ph, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), phlist ) )
  result_list = filter(x->x[1]!=nothing, pmap( ph->ph_evolve( p, funcs, ph, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), phlist ) )
  #result_list = filter(x->x[1]!=nothing, map( ph->ph_evolve( p, funcs, ph, max_tries, max_steps, use_lincircuit=use_lincircuit, use_mut_evolve=use_mut_evolve, print_steps=print_steps ), phlist ) )
  if length(result_list) == 0
    println("no results")
    return nothing
  end
  println("length(result_list): ",length(result_list))
  fail_fract = (length(phlist)-length(result_list))/length(phlist)
  steps_list = map(x->x[2], result_list)
  mean_steps = mean(steps_list)
  median_steps = median(steps_list)
  std_steps = std(steps_list)
  numactive_list = map(x->x[3], result_list )
  robust_list = map(x->x[4], result_list )
  mean_numactive = mean(numactive_list)
  mean_robust = mean(robust_list)
  total_failures_list = map(x->x[5], result_list )
  #println("total_failures_list: ",total_failures_list)
  first_fail = sum( total_failures_list[i][1] for i = 1:length(total_failures_list ))  # the number of phenotypes where evolution failed on the first try
  subsequent_fail = sum( sum( total_failures_list[i][2:end] ) for i = 1:length(total_failures_list ))  # the number of failures on the second and later tries
  df_row = [ p.numinputs, p.numinteriors, p.numlevelsback, fail_fract, mean_steps, median_steps, std_steps, mean_numactive, mean_robust, first_fail, subsequent_fail, total_failures_list ]
  push!(df,df_row)
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# MyInt: ",MyInt)
      print_parameters(f,p,comment=true)
      println(f,"# length(phlist): ",length(phlist))
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Evolve a circuit that maps to a given phenotype (Goal)
# There is another method of this function that evolves multiple circuits to the given phenotype
# max_tries is the maximum number of attempts using neutral_evolution to find the circuit
# max_steps is the maximum number of steps during a run of neutral_evolution()
# Note that there is no num_circuits_per_goal argument which distinguishes this version from the version that evolves multiple circuits per phenotype
# if no genotype that maps to goal is found, returns (nothing,nothing)
function pheno_evolve( p::Parameters, funcs::Vector{Func}, goal::Goal, max_tries::Int64, max_steps::Int64; 
    use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false )
  default_funcs(p)
  #Random.seed!(2)
  #println("length(funcs): ",length(funcs))
  steps = 0   # establish scope
  total_steps = 0   # establish scope
  total_failures = zeros(Int64,max_tries)
  nc = nothing # establish scope
  for i = 1:max_tries
    if use_lincircuit
      c = rand_lcircuit( p, funcs )
    else
      c = random_chromosome( p, funcs )
      #print_circuit(c)
    end         
    if !use_mut_evolve
      (nc,steps) = neutral_evolution( c, funcs, goal, max_steps, print_steps=print_steps )
      #(nc,steps) = neutral_evol( c, funcs, goal, max_steps, print_steps=print_steps );  println("neutral_evol")
    else
      (nc,steps,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [goal], funcs, max_steps ) 
    end
    total_steps += steps
    if steps < max_steps
      break
    end
    total_failures[i] += 1
  end
  if steps == max_steps
    println("pheno_evolve failed to evolve a circuit to goal: ",goal )
    return (nothing,nothing,nothing)
  end
  (nc,total_steps,total_failures)
end

# Evolve num_circuits_per_goal circuits that map to a given phenotype (Goal) goal.
# There is another method of this function that evolves one circuit to the given phenotype
# max_tries is the maximum number of attempts using neutral_evolution to find all circuits
#   Thus, max_tries needs to be greater than num_circuits_per_goal
# max_steps is the maximum number of steps during a run of neutral_evolution()
# Note that there is a num_circuits_per_goal argument which distinguishes this version from the version that evolves one circuit per phenotype
# Returns a list of 2-tuples (circuit,steps).  If no circuits are found, returns empty list
# if no genotype (circuit) that maps to goal is found, returns (nothing,nothing)
function pheno_evolve( p::Parameters, funcs::Vector{Func}, goal::Goal, num_circuits_per_goal::Int64, max_tries::Int64, max_steps::Int64; 
     use_lincircuit::Bool=false, use_mut_evolve::Bool=false, print_steps::Bool=false )
  if max_tries < num_circuits_per_goal
    error("Warning: set max_tries to be greater than num_circuits_per_goal since max_tries is for all calls to neutral_evolution().")
  end
  default_funcs(p)
  #println("length(funcs): ",length(funcs))
  steps = 0   # establish scope
  steps_per_circuit = 0  # The total number of evolution steps to evolve the current circuit
  nc = nothing # establish scope
  circuits_steps_list = use_lincircuit ? Tuple{LinCircuit,Int64}[] : Tuple{Chromosome,Int64}[]
  num_circuits_found = 0
  total_failures = zeros(Int64,max_tries)  # total_failures is not part of the return value of this function
  tries = 1
  while tries < max_tries && num_circuits_found < num_circuits_per_goal 
    tries += 1
    c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
    #print_circuit(c)
    if !use_mut_evolve
      (nc,steps) = neutral_evolution( c, funcs, goal, max_steps, print_steps=print_steps )
    else
      (nc,steps,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, [goal], funcs, max_steps ) 
    end
    steps_per_circuit += steps
    if steps < max_steps
      #print_circuit(nc)
      push!(circuits_steps_list, (nc,steps_per_circuit))
      steps_per_circuit = 0   # reset for next circuit
      num_circuits_found += 1
    else
      total_failures[tries] += 1
    end
  end
  if tries == max_tries 
    if num_circuits_found == 0
      println("pheno_evolve failed to evolve any circuits that map to goal: ",goal )
      return circuit_steps_list
    elseif num_circuits_found < num_circuits_per_goal
      println("pheno_evolved ",num_circuits_found," circuits that map to goal: ",goal )
      println("WARNING: less than num_circuits_per_goal which is ",num_circuits_per_goal)
    end
  end
  circuits_steps_list
end

# Starting with a chromosome c, epochal evolve numcircuits chromosomes that map to the phenotypes in phlist
function from_evolve( c::Chromosome, funcs::Vector{Func}, numcircuits::Int64, phlist::GoalList, max_tries::Int64, max_steps::Int64;
    use_mut_evolve::Bool=false, print_steps::Bool=false )::Tuple{Goal,Chromosome,Int64}    # The 3rd element of the tuple is steps
  result = Tuple{Goal,Chromosome,Int64}[]
  steps = nothing  # establish scope
  p = c.params
  for i = 1:numcircuits
    for ph in phlist
      for ttry = 1:max_tries
        if !use_mut_evolve
          (nc,steps) = neutral_evolution( c, funcs, ph, max_steps, print_steps=print_steps )
        else
          (nc,steps,worse,same,better,output,matched_goals,matched_goals_list) = mut_evolve( c, ph, funcs, max_steps ) 
        end
        if steps < max_steps
          push!(result,(ph,nc,steps))
          break  # from for ttry loop
        else  
          println("epochal evolution failed for phenotype ",ph," on ttry ",ttry)
        end
      end
    end
  end
  result
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
      #continue  # continue with current circuit c unchanged
    else
      println("c_new code:  ",circuit_code(c_new))
      println("c_dest_code: ",circuit_code(c_dest))
      new_dist = circuit_distance( c_new, c_dest )
      println("dist:     ",dist,"  new_dist: ",new_dist)
      if new_dist > dist
        step += 1
        #continue  # continue with current circuit c unchanged
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

# The name of the function neutral_evolution() is a misnomer since it actually implements what we call epochal evolution.
# Hu and Banzhaf call this an adaptive walk.
# Evolves a Chromsome or LinCircuit that maps to g starting with chromosome c.
# max_steps is the maximum number of evolutionary steps.
# If evolution hasn't succeeeded in max_steps, return nothing.
# insert_gate_prob is the probability of inserting a gate on a mutation of a chromosome.
# delete_gate_prob is similar for deleting a gate.
# select_prob is the probabability of accepting the new chromosome when there is an reduction in distance
# Similar to mut_evolve except that this takes a single goal instead of a goal list as an argument.
function neutral_evolution( c::Circuit, funcs::Vector{Func}, g::Goal, max_steps::Integer; print_steps::Bool=false,
      select_prob::Float64=1.0, gate_list::Vector{Int64}=Int64[], insert_gate_prob::Float64=0.0, delete_gate_prob::Float64=0.0, 
      save_acomplexities::Bool=false ) # Save acomplexity, evolvability, robustness of phenos produced by mutate_all on every step
  #println("neutral evolution started for goal: ",g)
  #if typeof(c) == Chromosome
  #  println("neutral evolution: ")
  #  print_circuit(c)
  #end
  default_funcs(c.params)
  p = c.params
  LinCirc = typeof(c) == LinCircuit ? :true : :false
  #if LinCirc println("LGP") else println("CGP") end
  step_list = Int64[] # Only used if save_acomplexities==true
  status_list = String[]  # Only used if save_acomplexities==true
  acomplexity_list = Float64[]  # Only used if save_acomplexities==true
  evolvability_list = Int64[]  # Only used if save_acomplexities==true
  robust_list = Float64[]   # Only used if save_acomplexities==true
  intersect_list = Int64[]  # Only used if save_acomplexities==true 
  #println("LinCirc: ",LinCirc,"  Ones: ",@sprintf("0x%08x",Ones),"  CGP.Ones: ",@sprintf("0x%08x",CGP.Ones))
  #println("numgates: ",c.params.numinteriors)
  #println("select_prob: ",select_prob)
  step = 0
  ov = output_values( c )
  #println("ov: ",ov,"  g: ",g)
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  while step < max_steps && ov != g
    step += 1
    #println("A step: ",step," current_distance: ",current_distance)
    new_c = deepcopy(c)
    if typeof(c) == Chromosome
      if length(gate_list) == 0
        (new_c,active) = mutate_chromosome!( new_c, funcs, insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
      else (new_c,active) = mutate_chromosome_gates!( new_c, funcs, gate_list )
      end  
    elseif typeof(c) == LinCircuit
      #println("c: ",c)
      new_c = mutate_circuit!( new_c, funcs )
    end
    new_ov = output_values( new_c )
    new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    #=
    if step < 20
      #print("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur dis: ",current_distance,"  new_dis: ",new_distance,"  " )
      #print_circuit(new_c)
    end
    =#
    #print("new_c: ")
    #print_circuit(new_c)
    #print("    c: ")
    #print_circuit(c)
    if new_ov == ov 
      c = new_c
      if save_acomplexities
        save_acomplexity( new_c, c, new_ov[1], g, step, "neutral", step_list, status_list, acomplexity_list, evolvability_list, robust_list, intersect_list )
      end
      if print_steps && step < 500
        print("step: ",step," is pheno neutral.  new_ov: ",new_ov,"  new_distance: ",new_distance,"  ")
        print_circuit(new_c)
      end
    elseif new_distance == current_distance
      c = new_c
      if save_acomplexities
        save_acomplexity( new_c, c, new_ov[1], g, step, "neutral", step_list, status_list, acomplexity_list, evolvability_list, robust_list, intersect_list )
      end
      if print_steps && step < 500
        print("step: ",step," is fitness neutral.  new_ov: ",new_ov,"  new_distance: ",new_distance,"  ")
        print_circuit(new_c)
      end
    elseif new_distance < current_distance   # improvement
      rnd = rand()
      if print_steps
        if rnd <= select_prob
          println("step: ",step,"  rnd: ",rnd,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance, "  new_c accepted")
        else
          println("step: ",step,"  rnd: ",rnd,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance, "  new_c not accepted")
        end
        print_circuit(new_c)
      end
      if save_acomplexities
        save_acomplexity( new_c, c, new_ov[1], g, step, "improve", step_list, status_list, acomplexity_list, evolvability_list, robust_list, intersect_list )
      end
      if rnd <= select_prob 
        c = new_c
        ov = new_ov
        current_distance = new_distance
        #println("B step: ",step," current_distance: ",current_distance)
      end
    else
      if print_steps && step <= 20
        print("step: ",step,"  new_output: ",new_ov," current distance: ",current_distance," new: ",new_distance,"  ")
        print_circuit(new_c)
      end 
    end
  end
  if save_acomplexities
    df = DataFrame( :step=>step_list, :status=>status_list, :acomplexity=>acomplexity_list, :evolvability=>evolvability_list, :robustness=>robust_list,
        :count_goals=>intersect_list )
    return df
  end
  if step == max_steps
    println("neutral evolution failed with ",step," steps for goal: ",g)
    return (c, step)
  else
    println("neutral evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end

# Example to compute evodict:    
# gdf = read_dataframe("../data/3_14_23/phnet_matrix3_14_23G.csv")   # for Parameters(3,1,8,4)
# edict = Dict{MyInt,Int64}()   # Dict{UInt16, Int64}()
# for i = MyInt(0):0x00ff edict[i] = gdf.d_evolvability[i+1] end
# Does a sequence of epochal evolutions (calls of neutral_evolution()) with each evolution starting from the result of the previous evolution
function neutral_evolution_glist( c::Circuit, funcs::Vector{Func}, gl::GoalList, max_steps::Integer, evodict::Dict{MyInt,Int64}=Dict{MyInt,Int64}(); 
      print_steps::Bool=false, select_prob::Float64=1.0 )
  P = c.params
  kdict = kolmogorov_complexity_dict( P, funcs )
  #println("neutral evolution started for goal: ",g)
  total_steps = 0
  failures = 0
  Tcomplexity_list = Float64[]
  p_evol_list = Float64[]
  cur_ch = deepcopy(c)
  for g in gl
    ( cur_ch, step ) = neutral_evolution( cur_ch, funcs, g, max_steps )
    total_steps += step
    result = (step == max_steps) ? "failed" : "succeeded"
    failures += (step == max_steps) ? 1 : 0
    Tcomplexity = complexity5(cur_ch)
    pheno = output_values( cur_ch )
    push!(Tcomplexity_list,Tcomplexity)
    p_evol = get( evodict, pheno[1], 0 )
    push!(p_evol_list,p_evol)
    println( "evolution $(result) for goal: ", g, "  Tononi complexity: ",@sprintf("%4.2f",Tcomplexity),"  K complexity: ", kdict[ output_values(cur_ch)[1] ],"  p_evolvability: ", p_evol )
  end
  ( total_steps, failures, Tcomplexity_list, p_evol_list )
end

function run_neutral_evolution_glist( c::Circuit, funcs::Vector{Func}, num_goals::Int64, Kcomp::Int64, hdist::Int64, max_steps::Integer, nreps::Int64, evodict::Dict{MyInt,Int64}=Dict{MyInt,Int64}(); 
      print_steps::Bool=false, select_prob::Float64=1.0 )
  mean_total_steps = 0.0
  mean_failures = 0.0
  mean_Tcomplexity_list = zeros(Float64, num_goals )
  mean_p_evol_list = zeros(Float64, num_goals )
  for i = 1:nreps
    plist = generate_phenotypes( P, funcs, Kcomp, hdist, num_goals )
    (total_steps, failures, Tcomplexity_list, p_evol_list ) = neutral_evolution_glist( random_chromosome(P,funcs), funcs, map(x->[x],plist), max_steps, evodict )
    mean_total_steps += total_steps
    mean_failures += failures
    mean_Tcomplexity_list .+= Tcomplexity_list
    mean_p_evol_list .+= p_evol_list
  end
  mean_total_steps /= nreps
  mean_failures /= nreps
  mean_Tcomplexity_list ./= nreps
  mean_p_evol_list ./= nreps
  ( mean_total_steps, mean_failures, mean_Tcomplexity_list, mean_p_evol_list )
end

# Returns a MyInt with nbits set to 1
# The bits set to 1 are in the sub-unsigned-integer corresponding to P.numinputs
# For example, if MyInt == UInt16 and P.numinputs == 3, then only bit chosen from bits 0 to 7 can be set to 1.
function rand_bit_word( P::Parameters, nbits::Int64=1 )
  total_nbits = 2^P.numinputs
  result = MyInt(0)
  for i = 1:nbits
    shift = rand(0:(total_nbits-1))
    #println("i: ",i,"  shift: ",shift)
    rand_bit_word = MyInt(1) << shift
    result = xor( result, rand_bit_word )
  end
  result
end

#function fitfunct_from_fitness_vector( P::Parameters, fitness::Vector{Float64} )::Function
#  fitfunct(x::MyInt) = fitness[x+1]
#end 

function testf( c::Circuit, fitfunct::Union{Vector{Float64},Function} )
  P = c.params
  println("Parameters: ",P)
  if typeof(fitfunct) <: Vector
    fitfunct = fitfunct_from_fitness_vector( P, fitfunct )
  end
  println("fitfunct(0x0004) ",fitfunct(0x0004))
end
   
# Uses epochal evolution based on fitness function represented as a vector to attempt to find a fitness non-decreasing path from cicuit c to a genotype of goal.
# Thus a requirement is that the fitness of circuit c is less than or equal to the fitness of goal.
### fitness is a vector of fitnesses indexed over all phenotypes (which limits the number of inputs to 4).
# fitfunct is either a vector of fitness indexed over MyInts, or a function from MyInts to Float64s.  MyInts describe single-output phenotypes
using GLM
function epochal_evolution_fitness( c::Circuit, funcs::Vector{Func}, g::Goal, fitfunct::Union{Vector{Float64},Function}, max_steps::Integer; 
      test_simplicity::Bool=false, print_steps::Bool=false )::Tuple{Union{Circuit,Nothing},Int64,Float64,Float64,Float64}
  function fitfunct_from_fitness_vector( P::Parameters, fitness::Vector{Float64} )::Function
    fitfunct(x::MyInt) = fitness[x+1]
  end 
  P = c.params
  if typeof(fitfunct) <: Vector
    fitfunct = fitfunct_from_fitness_vector( P, fitfunct ) 
  end
  src_ph = output_values( c )
  src_fit = fitfunct( src_ph[1] )
  dest_fit = fitfunct( g[1] )
  #println("src_ph: ",src_ph,"  src_fit: ",src_fit,"  dest_fit: ",dest_fit)
  @assert src_fit <= dest_fit
  kdict = test_simplicity ? kolmogorov_compexity_dict(P,funcs) : nothing
  local intercept_slope
  local predicted_Kcomp
  if test_simplicity
    pheno_range = MyInt(0):MyInt(2^2^P.numinputs-1)
    df = DataFrame(:goal=>map(x->[x],pheno_range),:Kcomp=>map(x->kdict[x],pheno_range), :half_unitation=>map(x->half_unitation(x,P),pheno_range))
    reg = GLM.lm(GLM.@formula(Kcomp ~ half_unitation), df)
    intercept_slope = GLM.coef(reg)
    predicted_Kcomp(half_unitation::Int64) = intercept_slope[1] + intercept_slope[2]*half_unitation  # Local function
    Kcomp_deviations = Float64[]
  end
  step = 0
  current_fitness = src_fit
  while step < max_steps && current_fitness < dest_fit
    step += 1
    #println("A step: ",step," current_fitness: ",current_fitness)
    new_c = deepcopy(c)
    if typeof(c) == Chromosome
      (new_c,active) = mutate_chromosome!( new_c, funcs )
    elseif typeof(c) == LinCircuit
      new_c = mutate_circuit!( new_c, funcs )
    end
    #println("output_values(new_c): ",output_values(new_c))
    #println("new_c params: ",new_c.params)
    new_fitness = fitfunct( output_values( new_c )[1] )
    #println("new_fitness: ",new_fitness)
    if new_fitness > current_fitness && new_fitness <= dest_fit  # improvement
      print_steps ? println("step: ",step,"  fitness improved from ",current_fitness," to ",new_fitness) : nothing
      c = new_c
      current_fitness = new_fitness
    elseif new_fitness > current_fitness && new_fitness > dest_fit  # improvement but higher fitness than goal
      print_steps ? println("step: ",step,"  fitness improved from ",current_fitness," to ",new_fitness," but fitness greater than goal fitness") : nothing
    elseif new_fitness == current_fitness  # neutral
      print_steps ? println("step: ",step,"  fitness neutral ",current_fitness) : nothing
      c = new_c
    else # fitness decrease
      print_steps ? println("step: ",step,"  fitness decrease from ",current_fitness," to ",new_fitness) : nothing
    end
    if test_simplicity
      ov = output_values(c)[1]
      push!(Kcomp_deviations, kdict[ov] - predicted_Kcomp( half_unitation( ov, P ) ) )
      #println("step: ",step,"  Kcomp_deviations[end]: ",Kcomp_deviations[end])
    end
  end # while
  if step == max_steps
    println("epochal evolution fitness failed with ",step," steps for goal: ",g)
    return (nothing, step, src_fit, dest_fit, 0.0 )
  else
    println("epochal evolution fitness succeeded at step ",step," for goal: ",g)
    @assert fitfunct( output_values( c )[1] ) == fitfunct( g[1] )
    if test_simplicity
      return (c, step, src_fit, dest_fit, mean(Kcomp_deviations) )
    else
      return (c, step, src_fit, dest_fit, 0.0 )
    end
  end
end

# Example run: run_epochal_evolution_fitness(P3, funcs, fitness, 1000, 2 )
function run_epochal_evolution_fitness( P::Parameters, funcs::Vector{Func}, fitfunct::Union{Vector{Float64},Function}, max_steps::Integer, nreps::Int64; 
    csvfile::String="", print_steps::Bool=false )::DataFrame
  if typeof(fitfunct) <: Vector
    fitfunct = fitfunct_from_fitness_vector( P, fitfunct )
  end
  df = DataFrame()
  df.numinputs = Int64[]
  df.numgates = Int64[]
  df.lb = Int64[]
  df.max_steps = Int64[]
  df.nreps = Int64[]
  df.src_fit = Float64[]
  df.dest_fit = Float64[]
  df.fit_diff = Float64[]
  df.dest_freq = Int64[]
  df.dest_Kcomp = Int64[]
  df.steps = Int64[]
  df.failures = Int64[]
  rdict = redundancy_dict(P,funcs)
  kdict = kolmogorov_complexity_dict(P,funcs)
  for i = 1:nreps
    (src,destgoal) = assign_src_dest( P, funcs, fitfunct )
    println("destgoal: ",destgoal)
    src_fit = fitfunct( output_values(src)[1] )
    dest_Kcomp = kdict[destgoal[1]]
    dest_freq = rdict[destgoal[1]]
    dest_fit = fitfunct( destgoal[1] )
    #println("src_fit: ",src_fit,"  dest_fit: ",dest_fit)
    (dch,steps) = epochal_evolution_fitness(src, funcs, destgoal, fitfunct, max_steps, print_steps=false )
    failure = (dch == nothing) ? 1 : 0
    df_row = (P.numinputs,P.numinteriors,P.numlevelsback,max_steps,nreps,src_fit,dest_fit,dest_fit-src_fit,dest_freq,dest_Kcomp,steps,failure)
    push!( df, df_row )
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs)
      print_parameters(f,P,comment=true)
      println(f,"# nreps: ",nreps)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Returns source chromosome and destination goal satisfying the condition fitness( src ) < fitness( destgoal )
function assign_src_dest( P::Parameters, funcs::Vector{Func}, fitfunct::Union{Vector{Float64},Function} )
  done = false
  src_fit = 1.0
  dest_fit = 0.0
  destgoal = nothing
  src = nothing
  done = false
  while !done
    src = random_chromosome(P,funcs)
    #src_fit = fitness[ output_values(src)[1]+1 ]
    src_fit = fitfunct( output_values(src)[1] )
    destgoal = randgoal( P )
    #dest_fit = fitness[ destgoal[1]+1 ]
    dest_fit = fitfunct( destgoal[1] )
    done = src_fit < dest_fit
  end
  #println("assign_src_dest: src_fit: ",src_fit,"  dest_fit: ",dest_fit)
  (src,destgoal)
end

function save_acomplexity( new_c::Circuit, prev_c::Circuit, ov::MyInt, g::Goal, step::Int64, status::String, step_list::Vector{Int64}, 
      status_list::Vector{String}, acomplexity_list::Vector{Float64}, evolvability_list::Vector{Int64}, robust_list::Vector{Float64}, intersect_list::Vector{Int64} )
  p = new_c.params
  funcs = default_funcs(p.numinputs)
  push!( step_list, step )
  push!( status_list, status )
  phenos = map( x->x[1], mutate_all( new_c, funcs, output_outputs=true ) )
  #println("phenos: ",phenos)
  #println("ph int g: ",length(findall( x->x==g[1], phenos )))
  push!( acomplexity_list, adami_complexity( phenos, p ) )
  #push!( acomplexity_list, mutual_information( phenos, fill( g[1], length(phenos ) ) ) ) 
  push!( evolvability_list, length(unique(phenos) ) )
  push!( robust_list, length( findall( x->x==ov, phenos ) )/length(phenos) )
  push!( intersect_list, length(findall( x->x==g[1], phenos )) )
end

# Computes phenotype evolvability, genotype evolvability, robustness, complexity, steps
# num_circuits is the number of circuits used to compute properties.
function geno_circuits( g::Goal, p::Parameters, funcs::Vector{Func}, num_circuits::Integer, max_steps::Integer, max_attempts::Integer )
  #funcs = default_funcs( p.numinputs )
  c = random_chromosome( p, funcs )
  sum_steps = 0
  circuit_list = Chromosome[]
  n_circuits = 0
  attempt = 0
  nc = nothing
  while attempt < max_attempts && n_circuits < num_circuits
    attempt += 1
    (nc,steps) = neutral_evolution( c, funcs, g, max_steps )
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

#=
# Computes phenotype evolvability, genotype evolvability, robustness, complexity, steps
# num_circuits is the number of circuits used to compute properties.
function geno_circuits( g::Goal, p::Parameters, funcs::Vector{Func}, num_circuits::Integer, max_steps::Integer, max_attempts::Integer )
  #funcs = default_funcs( p.numinputs )
  c = random_chromosome( p, funcs )
  sum_steps = 0
  circuit_list = Chromosome[]
  n_circuits = 0
  attempt = 0
  nc = nothing
  while attempt < max_attempts && n_circuits < num_circuits
    attempt += 1
    (nc,steps) = neutral_evolution( c, funcs, g, max_steps )
    sum_steps += steps
    if steps < max_steps
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
=#
 
# Added funcs to parameters on 4/28/23
function geno_properties( cl_ss::Tuple{Vector{Chromosome}, Int64}, funcs::Vector{Func} ) # Tuple of circuit list and sum_steps.
  (cl, sum_steps) = cl_ss
  #funcs = default_funcs( cl[1].params.numinputs )
  g = output_values(cl[1])
  sum_robust = 0.0
  sum_complexity = 0.0
  sum_evolvability = 0.0  # genotype evolability
  genotypes = Goal[]
  genotype_set = Set(Goal[])
  n_circuits = length(cl)
  i = 1
  for c in cl
    if !(output_values(c) == g)
      println("i: ",i,"  g: ",g,"  output_values(c): ",output_values(c),"  @assert output_values(c) == g failed in geno_properties() ")
      print_circuit(c)
    end
    #@assert output_values(c) == g
    #sum_robust += mutational_robustness( c, funcs )
    sum_robust += robustness( c, funcs )
    sum_complexity += complexity5(c) 
    genotypes = Set(mutate_all(c,funcs,output_outputs=true,output_circuits=false))
    sum_evolvability += length(genotypes)-1
    #println("len(genotypes): ",length(genotypes),"  genotypes: ",genotypes)
    genotype_set = union( genotype_set, genotypes )
    i += 1
  end
  (g, sum_robust/n_circuits, sum_complexity/n_circuits, sum_evolvability/n_circuits, length(genotype_set), sum_steps/n_circuits) 
end

# Computes the properites robustness, complexity, genotype evolability, phenotype evolvability, evolutonary steps
#    for the goals of gl.  The properties are based on the evoution of num_circuits circuits that map to the goal.
# max_attempts is the maximum number of calls to neutral_evolution() that are made in the attempt to evolve num_circuits circuits
#    that evolve to a goal 
# max_steps is the maximum number of evolutionary steps done by neutral_evoluion() during 1 attempt to evolve a circuit
#    that maps to a goal.
# Added funcs to parameters on 4/28/23
function geno_list_properties( gl::GoalList, p::Parameters, funcs::Vector{Func}, num_circuits::Integer, max_steps::Integer, max_attempts::Integer;
    csvfile::String="" )
  println("funcs: ",funcs)
  df = DataFrame()
  df.goal = Goal[]
  df.robustness = Float64[]
  df.complexity = Float64[]
  df.g_evolvability = Float64[]
  df.p_evolvability = Float64[]
  df.steps = Float64[]
  gp_list = Vector{Tuple{Goal,Float64,Float64,Float64,Float64,Float64}}[]
  gp_list = pmap( g->geno_properties( geno_circuits( g, p, funcs, num_circuits, max_steps, max_attempts ), funcs ), gl )
  #gp_list = map( g->geno_properties( geno_circuits( g, p, funcs, num_circuits, max_steps, max_attempts ), funcs ), gl )
  for gp in gp_list
    push!(df, gp) 
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs)
      print_parameters(f,p,comment=true)
      println(f,"# num_circuits: ",num_circuits)
      println(f,"# max_steps: ",max_steps)
      println(f,"# max_attempts: ",max_attempts)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Does epochal evolution using a 1+lambda strategy where lambda is specified in c.params.
function lambda_evolution( c::Chromosome, g::Goal, maxsteps::Integer, mutrate::Float64 )
  p = c.params
  p.mutrate = mutrate
  funcs = default_funcs(p.numinputs)
  poisson_lambda = p.mutrate*p.numinteriors
  X = Poisson( poisson_lambda )
  step = 0
  ov = output_values( c)
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  #print("ov: ",ov,"  cur_dist: ",current_distance,"  "); print_circuit(c)
  while step < maxsteps && ov != g
    chrome_list = Chromosome[]
    dist_list = Float64[]
    for i = 1:p.lambda
      step += 1
      num_mutations =  rand(X)
      new_c = deepcopy(c) 
      for m = 1:num_mutations
        mutate_chromosome!( new_c, funcs )
      end
      new_ov = output_values( new_c)
      new_dist = hamming_distance( new_ov, g, c.params.numinputs )
      #print("new_ov: ",new_ov,"  new_dist: ",new_dist,"  "); print_circuit(new_c)
      push!(chrome_list,new_c)
      push!(dist_list,new_dist)
    end  
    (best_dist,ind) = findmin( dist_list )
    if best_dist <= current_distance
      c = chrome_list[ind]
      ov = output_values( c )
      current_distance = hamming_distance( ov, g, c.params.numinputs )
      #print("ov: ",ov,"  cur_dist: ",current_distance,"  "); print_circuit(c)
    end
  end # while
  if step == maxsteps
    println("lambda evolution failed with ",step," steps for goal: ",g)
    return (c, step)
  else
    println("lambda evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end
# A simple no-frills pop_evolve().

function random_population( p::Parameters, popsize::Int64, funcs::Vector{Func} )
  pop = [ random_chromosome(p,funcs) for i = 1: popsize ]
end

function simple_pop_evolve( pop::Vector{Chromosome}, gl::GoalList, ngens::Int64, mutrate::Float64=1.0;
    prdebug::Bool=false )
  df = DataFrame()
  df.gen = Int64[]
  df.fract_optimal = Float64[]
  df.avg_pop_entropy = Float64[]
  df.mutual_inf = Float64[]
  funcs = default_funcs( p.numinputs )
  target = gl[1][1]
  for gen = 1:ngens
    fitness_vector = rescale_fitnesses([ pop[i].fitness for i = 1:popsize ])
    prdebug ? println("fit_vect: ",fitness_vector) : nothing
    prdebug ? println("gen: ",gen,"  max fit: ",maximum(fitness_vector)) : nothing
    propsel!( pop, fitness_vector, maxfit=findmax(fitness_vector)[1] )
    prdebug ? println("after propsel gen: ",gen) : nothing
    for c in pop
      #c.fitness = hamming_distance( gl, output_values(c), p.numinputs )
      c.fitness = fitness_funct( p, c, gl )
      prdebug ? print_circuit(c,include_fitness=true,include_robustness=false,include_pheno=true) : nothing
    end       
    #println("fract opt: ",@sprintf("%3.2f",fract_optimal_chromes( pop )),"  ave_ent: ",@sprintf("%3.2f",avg_pop_entropy( pop )))
    push!(df,[gen,fract_optimal_chromes(pop),avg_pop_entropy(pop),mut_info(pop,target,funcs)])
    sav_pop = deepcopy( pop )  
    for i in 1:popsize
      c = pop[i]
      if rand() <= mutrate
        pop[i] = c = mutate_chromosome!(deepcopy(c),funcs)[1]   
      end
    end
  end
  #pop
  df
end

function print_pop( pop::Vector{Chromosome} )
  for c in pop
    print_circuit(c,include_fitness=true,include_robustness=false,include_pheno=true)
  end
end

function fract_optimal_chromes( pop::Vector{Chromosome} )
  count = 0
  for c in pop
    if c.fitness == 1.0
      count += 1
    end
  end
  count/length(pop)
end

# Copied from Pop_evolve.jl.  Remove
function rescale_fitnesses( fit_vect::Vector{Float64} )
  fit_min = minimum( fit_vect )
  #fit_min = quantile( fit_vect, 0.20 )
  fit_max = maximum( fit_vect )
  frange = fit_max - fit_min
  if frange > 0.0
    return [ (f-fit_min)/frange for f in fit_vect ]
  else
    return fit_vect
  end
end        

function fitness_funct( p::Parameters, c::Chromosome, gl::GoalList; mut_inf_matrix::Matrix{Float64}=zeros(Float64,1,1) )
  ov = output_values(c)
  #println("size(mut_inf_matrix)[1]: ",size(mut_inf_matrix)[1])
  if size(mut_inf_matrix)[1] <= 1
    maximum( 1.0-hamming_distance( g, ov, p.numinputs ) for g in gl )
  else  
    maximum( mut_inf_matrix[ov[1]+1,g[1]+1] for g in gl )
  end
end

function run_to_rand_phenos( param_list::Vector{Parameters}, funcs::Vector{Func}, ngoals::Int64, ngoalreps::Int64,  max_tries::Int64, max_steps::Int64;
    csvfile::String="" )
  phlist = randgoallist( ngoalreps*ngoals, paramlist[1], repetitions=ngoalreps )
  println("phlist: ",phlist)
  df = DataFrame( :numinputs=>Int64[], :numgates=>Int64[], :levsback=>Int64[], :steps=>Float64[], :tries=>Float64[], :Tcmplx=>Float64[] )
  for p in param_list
    default_funcs(p)
    ( steps, tries, Tcmplx ) = run_evolve_to_pheno( p, funcs, phlist, max_tries, max_steps )
    push!(df,(p.numinputs, p.numinteriors, p.numlevelsback, steps, tries, Tcmplx ))
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(param_list[1].numinputs))
      println(f,"# numinputs_list: ",map(i->param_list[i].numinputs,1:length(param_list)))
      println(f,"# numinteriors_list: ",map(i->param_list[i].numinteriors,1:length(param_list)))
      println(f,"# ngoals: ",ngoals)
      println(f,"# ngoalreps: ",ngoalreps)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

#=
function run_to_rand_phenos_mt( param_list::Vector{Parameters}, funcs::Vector{Func}, ngoals::Int64, ngoalreps::Int64,  max_tries::Int64, max_steps::Int64;
    csvfile::String="" )
  phlist = randgoallist( ngoalreps*ngoals, param_list[1], repetitions=ngoalreps )
  println("phlist: ",phlist)
=#

function run_to_rand_phenos_mt( param_list::Vector{Parameters}, funcs::Vector{Func}, ngoals::Int64, ngoalreps::Int64,  max_tries::Int64, max_steps::Int64;
    csvfile::String="" )
  phlist = randgoallist( ngoalreps*ngoals, param_list[1], repetitions=ngoalreps )
  println("phlist: ",phlist)
  df = DataFrame( :numinputs=>Int64[], :numgates=>Int64[], :levsback=>Int64[], :ngoals=>Int64[], :ngoalreps=>Int64[],  
      :steps=>Float64[], :tries=>Float64[], :Tcmplx=>Float64[], :robustness=>Float64[], :tries_list=>Vector{Int64}[] ) 
  for p in param_list
    default_funcs(p)  # sets the global variable Ones
    ( steps, tries, Tcmplx, robustness, tries_list ) = run_evolve_to_pheno_mt( p, funcs, phlist, max_tries, max_steps )
    push!(df,(p.numinputs, p.numinteriors, p.numlevelsback, ngoals, ngoalreps, steps, tries, Tcmplx, robustness, tries_list ))
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(param_list[1].numinputs))
      println(f,"# numinputs_list: ",map(i->param_list[i].numinputs,1:length(param_list)))
      println(f,"# numinteriors_list: ",map(i->param_list[i].numinteriors,1:length(param_list)))
      println(f,"# ngoals: ",ngoals)
      println(f,"# ngoalreps: ",ngoalreps)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function run_evolve_to_pheno_mt( p::Parameters, funcs::Vector{Func}, phlist::GoalList, max_tries::Int64, max_steps::Int64 )
  println("run_evolve_to_pheno_mt() started")
  nphenos = length(phlist)
  println("nphenos: ",nphenos)
  robust = 0.0
  #ttry_list = [ Atomic{Float64}(0.0) for i= 1:nphenos]
  #steps_list = [ Atomic{Float64}(0.0) for i= 1:nphenos]
  #Tcomplexity_list = [ Atomic{Float64}(0.0) for i= 1:nphenos]
  ttry_list = [ Atomic{Int64}(0) for i= 1:max_tries]
  ttries = Atomic{Int64}(0.0)
  ssteps = Atomic{Int64}(0.0)
  Tcomplexity = Atomic{Float64}(0.0)
  zrobustness = Atomic{Float64}(0.0)
  Threads.@threads for i = 1:nphenos
    ( nc, steps, tries ) = evolve_to_pheno( p, funcs, phlist[i], max_tries, max_steps )
    Threads.atomic_add!( ssteps, steps )
    Threads.atomic_add!( ttries, tries )
    Threads.atomic_add!( ttry_list[tries], 1 )
    Tcmplx = (p.numinteriors <= 20) ? complexity5(nc) : 0.0
    Threads.atomic_add!( Tcomplexity, Tcmplx )
    #println("before robust computation")
    robust =robustness(nc,funcs)
    Threads.atomic_add!( zrobustness, robust )
    #println("robust: ",robust)
  end  
  println("run_evolve_to_pheno_mt() finished")
  ( ssteps[]/nphenos, ttries[]/nphenos, Tcomplexity[]/nphenos, zrobustness[]/nphenos, map(tt->tt[],ttry_list) )
end

function run_evolve_to_pheno( p::Parameters, funcs::Vector{Func}, phlist::GoalList, max_tries::Int64, max_steps::Int64 )
  ttry_list = Int64[]
  steps_list = Int64[]
  Tcomplexity_list = Float64[]
  robustness_list = Float[]
  result_list = pmap( ph->evolve_to_pheno( p, funcs, ph, max_tries, max_steps ), phlist )
  #result_list = map( ph->evolve_to_pheno( p, funcs, ph, max_tries, max_steps ), phlist )
    #(nc, ttry, steps) = evolve_to_pheno( p, funcs, ph, max_tries, max_steps )
  for result in result_list
    push!(steps_list,result[2])
    push!(ttry_list,result[3])
    if p.numinteriors <= 20
      push!(Tcomplexity_list,complexity5(result[1]))
    else
      push!(Tcomplexity_list,0.0)
    end
  end
  (mean(steps_list), mean(ttry_list), mean(Tcomplexity_list))
end    

# Duplicates functionality of function pheno_evolve()
function evolve_to_pheno( p::Parameters, funcs::Vector{Func}, ph::Goal, max_tries::Int64, max_steps::Int64; c::Chromosome=random_chromosome(Parameters(1,1,0,1),default_funcs(p)) )
  if c.params.numinputs == 1
    c = random_chromosome(p,funcs)
  end
  ttry = 1
  (nc,steps) = neutral_evolution( c, funcs, ph, max_steps)
  total_steps = steps  # Never used
  while ttry < max_tries && steps == max_steps
    (nc,steps) = neutral_evolution( c, funcs, ph, max_steps )
    total_steps += steps
    ttry += 1
  end
  if ttry == max_tries && steps == max_steps
    return nothing
  end
  return (nc,steps,ttry)  # steps is the number of steps of the last successful run
end

function k_complexity_adaptive_walk( p::Parameters, funcs::Vector{Func}, g::Goal, max_ev_steps::Int64, max_tries::Int64, walk_steps::Int64 )
  c = random_chromosome( p, funcs )
  nc = neutral_evolution( c, funcs, g, max_ev_steps )
end
  
function sampling_evolution( c::Circuit, funcs::Vector{Func}, g::Goal, cdf::DataFrame, max_steps::Integer ) 
  @assert typeof(c) == Chromosome
  @assert "circuit_ints" in names(cdf)
  p = c.params
  circ_int_counter = zeros( Int64, 2^2^p.numinputs )
  ph = output_values(c)  # current phenotype
  step = 0
  while step < max_steps && ph !== goal
    phi = Int64(ph[1])+1
    len = length(cdf.circuit_ints[phi])
    j = circ_int_counter[phi] % len
    circ_int_counter[phi] += 1
    new_c = int_to_chromosome( cdf.circuit_ints[circ_int_counter[phi]][j+1] )
    ph = output_values(c)  # current phenotype
  end
end

# Does neutral evolution starting with Chromosome trial and goal Chromosome target.
# Each step uses mutate_all() to look for the next Chromosome.
# If no chromosome closer to target is found, the evolution stops.
# Returns a list of triples where each triple is a (phenotype, chromosome, distance )
# Seems to always reduce distance by 1.
function directed_neutral_evolution( trial::Chromosome, funcs::Vector{Func}, target::Chromosome, max_steps::Int64 )
  p = trial.params
  triples_list = nothing
  fmin = nothing
  initial_distance = geno_distance( trial, target, funcs )
  #println("initial distance: ",initial_distance)
  current_distance = initial_distance
  ch = deepcopy(trial)
  for step = 1:max_steps
    ch_ph = output_values(ch)
    #println("step: ",step,"  ch_ph: ",ch_ph)
    (outputs_list, ch_list) = mutate_all( ch, funcs, output_outputs=true, output_circuits=true )
    pairs_list = filter(x->x[1][1]==ch_ph[1], [ (outputs_list[i],ch_list[i]) for i = 1:length(ch_list) ])
    triples_list = [ (pairs_list[i][1], pairs_list[i][2], geno_distance( pairs_list[i][2], target, funcs )) for i = 1:length(pairs_list) ]
    distances = map(x->x[3],triples_list)
    new_distance = minimum(distances)
    fmin = (new_distance, rand(findall(x->x==new_distance,distances)))
    #println("fmin: ",fmin)
    new_ch = pairs_list[fmin[2]][2]
    #println("new_distance: ",new_distance,"  new_ch: ",new_ch)
    if new_distance >= current_distance
      current_distance = new_distance
      break
    else
      current_distance = new_distance
      ch = deepcopy(new_ch)
    end
  end
  #fmin
  findmin(map(x->x[3],triples_list))[1]
end
