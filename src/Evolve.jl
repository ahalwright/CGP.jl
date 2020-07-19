export components_matched, goals_matched_exact, goals_matched_hamming, next_chromosome!, mut_evolve, mut_reduce_numactive
export random_neutral_walk, match_score, findmaxall, findminall, findmaxall_1, findminall_1 


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
    mxscores = [ (sum( H[i,p[i]] for i = 1:nc ),[p[i] for i = 1:nc]) for p in P ]
    findmaxall(mxscores)[1][1]
  else # number of exact matches plus Hamming score for best partial match
    mxscores = [ maxscore([H[i,p[i]] for i = 1:nc]) for p in P ]
    findmaxall(mxscores)[1]
  end
end
#=
function match_score( output::Vector{MyInt}, goal::Vector{MyInt},numinputs::Int64)
  nc = length(output)  # number of components
  #println("output: ",output,"  goal: ",goal)
  @assert nc == length(goal)
  H = [ hamming_distance(output[i],goal[j],numinputs) for i = 1:nc, j=1:nc]  
  #println("H: ",H)
  P = permutations(collect(1:nc))
  mxscores = [ maxscore([H[i,p[i]] for i = 1:nc]) for p in P ]
  #println("mxscores: ",mxscores)
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
function goals_matched_exact( output::Vector{MyInt}, goallist::GoalList, numinputs::Int64 )
  #println("exact")
  matched_goals = map( g->components_matched(output,g), goallist )
  #println("matched_goals: ",matched_goals)
  fm = findmaxall( map(length, matched_goals ))
  #println("fm: ",fm)
  (Float64(fm[1]), fm[2], matched_goals[fm[2]])
end

# Note that fitness does not take robustness into account
function goals_matched_hamming( output::Vector{MyInt}, goallist::GoalList, numinputs::Int64; avgfitness::Bool=false  )
  #println("hamming")
  goallist_scores = map( g->match_score(output,g,numinputs,avgfitness=avgfitness), goallist )
  matched_goals = map( g->components_matched(output,g), goallist )
  #println("goallist_scores: ",goallist_scores)
  #println("matched_goals: ",matched_goals)
  (best_score,best_score_list) = findmaxall( goallist_scores )
  (best_score, best_score_list, matched_goals[best_score_list])
end

# Removed robust_sel and active_only keyword args on 6/6/20
function next_chromosome!(c::Chromosome, goallist::GoalList, funcs::Vector{Func}, fitness::Float64=0.0; 
      use_robustness::Bool=false, hamming_sel::Bool=true, active_only::Bool=false, num_mutations::Int64=1,
      avgfitness::Bool=false )
  context = construct_context(c.params.numinputs)
  new_c = deepcopy(c)
  for _ = 1:num_mutations
    mutate_chromosome!( new_c, funcs )
  end
  output = execute_chromosome(new_c,context)
  #println("output: ",output)
  goals_matched = hamming_sel ? goals_matched_hamming : goals_matched_exact
  ( new_fitness, matched_goals, matched_goals_list ) = goals_matched( output, goallist, c.params.numinputs,
      avgfitness=avgfitness )
  #println("new fitness: ",new_fitness)
  new_c = new_fitness >= fitness ? new_c : c
  if use_robustness
    mut_robust = mutational_robustness( new_c, funcs, active_only=active_only )
    new_c.robustness = mut_robust
    #println("new_fitness: ",new_fitness,"  n_mut_robust: ",mut_robust)
  end
  new_c.fitness = new_fitness
  return ( new_c, matched_goals, matched_goals_list )
end

# Does single-individual neutral evolution.  
# Does subgoal evolution where the number of matched components of a goal is maximized if exact=true
# Does Hamming evolution where the Hamming distance of matched components of a goal is minimized if exact=false
# Removed robust_sel and active_only keyword args on 6/6/20
function mut_evolve( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, max_steps::Integer;
      hamming_sel::Bool=true, use_robustness::Bool=false, num_mutations::Int64=1, print_improvements::Bool=false,
      avgfitness::Bool=false )
  #println("mut evolve avgfitness: ",avgfitness)
  #print_chromosome(c)
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
  (c, matched_goals, matched_goals_list ) = next_chromosome!(c, goallist, funcs, 
      hamming_sel=hamming_sel, use_robustness=use_robustness, num_mutations=num_mutations, avgfitness=avgfitness ) 
  output = output_values(c)   # Executes c if it has not already been executed
  #println("output: ",output)
  #println("trunc(c.fitness): ",trunc(c.fitness))
  #while step < max_steps && new_fitness < c.params.numoutputs
  while step < max_steps && trunc(c.fitness) < c.params.numoutputs
    #println("step: ",step,"  output: ",output,"   ")
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
    if step < max_steps
      (c, matched_goals, matched_goals_list ) = 
            next_chromosome!(c, goallist, funcs, hamming_sel=hamming_sel, use_robustness=use_robustness, 
            num_mutations=num_mutations, avgfitness=avgfitness ) 
      #=
      (c, matched_goals, matched_goals_list ) = 
            next_chromosome!(c, goallist, funcs, fitness, hamming_sel=hamming_sel, 
            use_robustness=use_robustness, avgfitness=avgfitness ) 
      =#
      output = output_values(c)   # Executes c if it has not already been executed
    end
  end
  if step == max_steps
    println("mut_evolve finished at step limit ",max_steps," with fitness: ", c.fitness ) 
  else
    #println("matched_goals: ",matched_goals,"  matched_goals_list: ",matched_goals_list)
    println("mut_evolve finished in ",step," steps with fitness: ", c.fitness )
  end
  if orig_c.fitness > c.fitness
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
    
    try
      @assert sort(output) == sort(goallist[matched_goals[1]])
    catch
      println("sort(output): ",sort(output),"  sort(goallist[matched_goals[1]]): ",sort(goallist[matched_goals[1]]))
    end
    
  end
  (c,step,worse,same,better,output,matched_goals,matched_goals_list)
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
    
# Change the fields of chromosome c to be the fields of chromosom c_to_copy
function copy_chromosome!( c::Chromosome, c_to_copy::Chromosome )
  c.params = c_to_copy.params
  c.inputs = c_to_copy.inputs
  c.interiors = c_to_copy.interiors
  c.outputs = c_to_copy.outputs
  c.interiors = c_to_copy.interiors
end
#=
mutable struct Int_node
  func::Func
  inputs::Vector{Int64}
end

mutable struct Input_node
  index::Int64
end

mutable struct Output_node
  input::Integer
end

# Example call: build_chromosome( [Input_node(1),Input_node(2)], [Int_node(OR,[1,2]),Int_node(AND,[2,3])],[Output_node(4)])
function build_chromosome( input_nodes::Vector{Input_node}, interior_nodes::Vector{Int_node}, output_nodes::Vector{Output_node} )
  num_in = length(input_nodes)
  num_ints = length(interior_nodes)
  num_outs = length(output_nodes)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=num_ints+num_outs )
  in_nodes = [InputNode(in_node.index) for in_node in input_nodes]
  int_nodes = [InteriorNode(int_node.func, int_node.inputs) for int_node in interior_nodes]
  out_nodes = [OutputNode(outnode.input) for outnode in output_nodes]
  Chromosome( p, in_nodes, int_nodes, out_nodes, 0.0, 0.0 )
end
=#

# Example call:  build_chromosome((1,2), ((OR,[1,2]),(AND,[2,3])),(4,))
function build_chromosome( inputs::Tuple, ints::Tuple, outs::Tuple )
  num_in = length(inputs)
  num_ints = length(ints)
  num_outs = length(outs)
  p = Parameters( numinputs=num_in, numoutputs=num_outs, numinteriors=num_ints, numlevelsback=num_ints+num_outs )
  in_nodes = [InputNode(in_index) for in_index in inputs]
  int_nodes = [InteriorNode(int_pair[1], int_pair[2]) for int_pair in ints]
  out_nodes = [OutputNode(out_index) for out_index in outs]
  Chromosome( p, in_nodes, int_nodes, out_nodes, 0.0, 0.0 )
end

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
