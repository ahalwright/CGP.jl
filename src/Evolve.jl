# Functions to implement mutational and combinational evolution

export mut_evolve_exact, mut_evolve_subgoal

#=
goal1 = MyInt[ 0xb, 0xe ]
goal2 = MyInt[ 0x6, 0x9 ]
goal3 = MyInt[ 0xa, 0xe ]
goallist2 = [goal1, goal2,goal3]
=#

# Do mutational evolution starting with chromosome c
# For each of the up to max_steps iterations of the loop, chromosome c is mutated to new_c, and
#    one of the following alternatives happens:
#  1.  new_c computes a new logic function which is not a goal
#      In this case, new_c is discarded
#  2.  new_c computes the same logic function (which is not a goal)
#      In this case, c becomes new_c, and a message is printed
#  3.  new_c exactly computes a goal of goallist
#      In this case, the iteration stops and the function returns a value described below.
# If a chromosome is found that exactly computes a goal in goallist, returns a triple of:
#  1.  The goal from goallist that is computed
#  2.  The previous logic function
#  3.  The original logic function computed by c
# Comment:  Should return new_c, the chromosome that computes the goal
function mut_evolve_exact( c::Chromosome, context::Vector{MyInt}, goallist::GoalList, max_steps::Integer )
  orig_c = deepcopy(c)
  orig_result = execute_chromosome( c, context )
  prev_c = deepcopy(c)
  prev_result = orig_result
  new_c = mutate_chromosome!( c, default_funcs() )
  new_result = execute_chromosome( new_c, context )
  if goal_check( new_result, goallist )
    println("goal found after 0 steps")
    return (new_result,prev_result,orig_result)
  end
  #println("orig_result: ",orig_result,"  new_result: ",new_result)
  step = 1
  #while !goal_check( execute_chromosome(new_c,context), tmp_goallist ) && step < max_steps
  while( (cc = chrome_check( new_result, prev_result, goallist )) != 2) && step < max_steps
    if cc == 1   # found a new chromosome with the same result
      prev_result = new_result
      println("mutated chromosome with same result found after ", step, " steps")
    end
    new_c = mutate_chromosome!( c, default_funcs() )
    new_result = execute_chromosome(new_c,context)
    #println("prev_result: ",prev_result,"  new_result: ",new_result)
    step += 1
  end
  if step == max_steps
    println("no goal found in ", max_steps, " steps" )
  else
    #cc = chrome_check( new_result, orig_result,  goallist ) 
    println("cc: ",cc)
    if cc == 2
      println("goal found after ", step, " steps")
    elseif cc == 1
      println("this should not happen")
      println("mutated chromosome with same result found after ", step, " steps")
    end
    #print_chromosome( new_c )
  end
  (new_result,prev_result,orig_result)
end

function chrome_check( new_result::Goal,previous_result::Goal, goallist::GoalList )
  if new_result == previous_result
    return 1
  elseif goal_check( new_result, goallist )
    return 2
  else 
    return 0
  end
end

#=
function mut_evolve_exact( c::Chromosome, context::Vector{MyInt}, goallist::GoalList, max_steps::Integer )
  new_c = c
  new_result = execute_chromosome( new_c, context )
  orig_c = deepcopy(c)
  orig_result = new_result
  prev_c = deepcopy(c)
  prev_result = orig_result
  #println("orig_result: ",orig_result,"  new_result: ",new_result)
  step = 1
  #while !goal_check( execute_chromosome(new_c,context), tmp_goallist ) && step < max_steps
  while( (cc = chrome_check( new_result, prev_result, goallist )) != 2) && step < max_steps
    if cc == 1   # found a new chromosome with the same result
      prev_result = new_result
      println("mutated chromosome with same result found after ", step, " steps")
    end
    new_c = mutate_chromosome!( c, default_funcs() )
    new_result = execute_chromosome(new_c,context)
    #println("step: ",step,"  prev_result: ",prev_result,"  new_result: ",new_result)
    step += 1
  end
  if step == max_steps
    println("no goal found in ", max_steps, " steps" )
  else
    #cc = chrome_check( new_result, orig_result,  goallist ) 
    println("cc: ",cc)
    if cc == 2
      println("goal found after ", step, " steps")
    elseif cc == 1
      println("this should not happen")
      println("mutated chromosome with same result found after ", step, " steps")
    end
    #print_chromosome( new_c )
  end
  (new_result,prev_result,orig_result)
end
=#
