
# Combines two 1-output circuits c1 and c2 with the same paramterers and number of gates into a 2-output circuit
# The output values of the result is the 2-output goal whose components are the output values of c1 and of c2.
function combine_circuits( c1::Chromosome, c2::Chromosome )
  # Two functions used only within this function
  inp1(i,numinputs) = i <= numinputs ? i : 2*(i-numinputs)+numinputs-1
  inp2(i,numinputs) = i <= numinputs ? i : 2*(i-numinputs)+numinputs
  @assert c1.params.numinputs == c2.params.numinputs
  @assert c1.params.numoutputs == 1
  @assert c2.params.numoutputs == 1
  if length(c1.interiors) < length(c2.interiors)
    println("  c1ints: ",length(c1.interiors),"  c2ints: ",length(c2.interiors))
    ctemp = c1
    c1 =c2
    c2 = ctemp
  end
  #println("c1ints: ",length(c1.interiors),"  c2ints: ",length(c2.interiors))
  p = Parameters( c1.params.numinputs, 2, c1.params.numinteriors+c2.params.numinteriors,
    c1.params.numlevelsback+c2.params.numlevelsback )
  inputs = c1.inputs
  interiors = InteriorNode[]
  new_numinteriors = 2*length(c2.interiors) 
  #println("new_numints: ",new_numinteriors)
  #print_parameters(p)
  for i = 1:length(c2.interiors)
    int_node1 = deepcopy(c1.interiors[i])
    input1 = inp1(c1.interiors[i].inputs[1],p.numinputs)
    input2 = inp1(c1.interiors[i].inputs[2],p.numinputs)
    int_node1.inputs = [input1,input2]
    #println("i: ",i,"  int_node1: ",int_node1)
    push!(interiors,int_node1)
    int_node2 = deepcopy(c2.interiors[i])
    input1 = inp2(c2.interiors[i].inputs[1],p.numinputs)
    input2 = inp2(c2.interiors[i].inputs[2],p.numinputs)
    int_node2.inputs = [input1,input2]
    #println("i: ",i,"  int_node2: ",int_node2)
    push!(interiors,int_node2)
  end
  #println("interiors: ",interiors)
  outputs = [OutputNode(p.numinputs+new_numinteriors-1),OutputNode(p.numinputs+new_numinteriors)]
  Chromosome( p, c1.inputs, interiors, outputs, 0.0, 0.0 )
end

# Starting with circuit c (and its parameters), try to evolve a circuit whose output is the same as c
#   and whose number of active nodes is reduced.
# Version of 5/7/21
function evolve_reduce_numactive( c::Chromosome, maxsteps::Int64; num_mutations::Int64=1, funcs::Vector{Func}=default_funcs(c.params.numinputs)  )
  p = c.params
  goal = output_values(c)
  current_num_active = number_active( c )
  for i = 1:maxsteps
    new_c = c
    for _ = 1:num_mutations
      new_c = deepcopy(c)
      mutate_chromosome!( new_c, funcs )
    end
    if output_values(new_c) == goal && number_active( new_c ) <= current_num_active
      c = new_c
      if number_active( c ) < current_num_active
        println("number active reduced from ",current_num_active," to ",number_active( new_c )," on step ",i,".")
        current_num_active = number_active( c ) 
      end
    end
  end
  (c, current_num_active )
end


# Use neutral evolution to find a chromosome c that maximizes funct(c)
# funct( c::Chromosome ) returns a Float64.
function evolve_function( funct, p::Parameters, funcs::Vector{Func}, max_steps::Int64;
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

# Attempts to reduce the number active of the chromosome c by mutational neutral evolution
# Older version written prior to 7/18/20.
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

