# Functions to implement goal functions and testing for goal functions

export Goal, GoalList, randgoal, randgoallist, rand_env_goallist, env_goal, goal_count, goal_check, my_test_goal

# Moved these to aliases.jl
#const Goal =  Vector{MyInt}
#const GoalList = Vector{Goal}

#=
goal1 = MyInt[ 0xb, 0xe ]
goal2 = MyInt[ 0x6, 0x9 ]
goal3 = MyInt[ 0x1, 0xe ]
goal4 = MyInt[ 0x1, 0xe ]
goal5 = MyInt[ 0x3, 0x9 ]
goallist2 = [goal1, goal2, goal3 ]

goallist3 = [0x50, 0x05, 0xf5, 0xa0, 0xe0, 0x0e, 0xf1, 0x80, 0x08, 0xfa]
=#

function component_type( numinputs::Int64 )
  if numinputs == 2
    comp_type = collect(0x00:0x0f)
  elseif numinputs == 3
    comp_type = collect(0x00:0xff)
  elseif numinputs == 4
    comp_type = UInt16
  elseif numinputs == 5
    comp_type = UInt32
  elseif numinputs == 6
    comp_type = UInt64
  elseif numinputs == 7
    comp_type = UInt128
  elseif numinputs > 7
    error("randgoal doesn't work for numinputs > 7")
  end
  comp_type
end

# generate a random goal assuming that MyInt is set correctly for numinputs
# Creates an inexacterror if conversion to MyInt doesn't work
function randgoal(numinputs::Int64, numoutputs::Int64 ) 
  map(x->convert( MyInt, x), rand( component_type( numinputs ), numoutputs ) )
end

# generate a random goallist of length ngoals
function randgoallist(ngoals::Int64, numinputs::Int64, numoutputs::Int64; repetitions::Int64=1)
  result = Vector{MyInt}[]
  if ngoals % repetitions != 0
    error("ngoals=",ngoals," must be a multiple of repetitions=",repetitions," in randgoallist")
  end
  #println("  collect(1:repetitions.ngoals): ",collect(1:repetitions:ngoals))
  for i = 1:repetitions:ngoals
    goal = rand( component_type(numinputs), numoutputs )
    for j = 1:repetitions
      #println("(i,j): ",(i,j)," i+j-1: ",i+j-1)
      #result[i+j-1] = goal
      push!(result,goal)
    end
  end
  result
end     

function rand_env_goallist( ngoals::Int64, numinputs::Int64, numoutputs::Int64, 
      repetitions::Int64, num_flipped_bits::Int64 )
  result = Vector{MyInt}[]
  if numoutputs % repetitions != 0
    error("numoutputs=",numoutputs," must be a multiple of repetitions=",repetitions," in rand_env_goallist")
  end
  f_bits_list = flipped_bits_list( num_flipped_bits, repetitions-1, numinputs )
  #println("f_bit_list: ",f_bits_list)
  for i = 1:ngoals
    #g = rand( component_type(numinputs), div( numoutputs, repetitions ) )
    g = randgoal( numinputs, div(numoutputs,repetitions) )
    goal = env_goal( g, f_bits_list )
    Base.push!(result,goal)
  end
  result
end 

function flipped_bits_list( num_flipped_bits::Int64, repetitions::Int64, numinputs::Int64 )
  one = MyInt(1)
  ones = one
  for i = 1:(num_flipped_bits-1)
    ones = ones<<1 | one
  end 
  fb = ones << (2^numinputs-num_flipped_bits)
  #Printf.@printf("flipped_bits_list: ones:  0x%04x  fb: 0x%04x\n",ones,fb)
  [ fb >> (num_flipped_bits*i) for i = 0:(repetitions-1) ]
end

function env_goal( goal::Goal, flipped_bits_list::Vector{MyInt} )
  num_repeats = length(flipped_bits_list)
  new_goal = Goal()
  for c in goal
    Base.push!( new_goal, c )
    for i = 1:num_repeats
      Base.push!( new_goal, apply_condition( c, flipped_bits_list[i] ))
    end
  end
  #println("new_goal: ",new_goal)
  new_goal
end
  
function apply_condition( component::MyInt, flipbits::MyInt )
  neg_flipbits = NOT.func( flipbits )
  neg_component = NOT.func( component )
  (neg_component & flipbits ) | (component & neg_flipbits )
end
  

# Returns a 2-tuple of:
#    1.  the  maximum number of matching subgoals of testgoal to components of goallist
#    2.  the  index of the goal of goallist that gave this maximum
# If testgoal equals all subgoals of a goal of goallist, then returns numouputs and the index of the goal in goallist
# If no component of testgoal matches the corresponding component of any goal in goallist, return 0
function goal_count( testgoal::Goal, goallist::GoalList )
   findmax( map( x -> my_test_goal( x, testgoal ), goallist ))[1] 
end

# Returns true if testgoal is equal to any of the goals in goallist defined above
function goal_check( testgoal::Goal, goallist::GoalList )
  length( filter( x -> testgoal==x, goallist )) > 0
end

# Returns the number of components of testGoal that are equal to the corresponding component of goal.
#function my_test_goal( testGoal::Goal, goal::Goal )
#  length( filter( x -> x == 0, testGoal .⊻ goal ))
#end
function my_test_goal( testGoal::Goal, goallist::GoalList )
  achieved = 0
  for g in goallist
    ga = length( filter( x -> x == 0, testGoal .⊻ g ))
    if ga > achieved
      achieved = ga
    end
    println("my_test_goal: achieved: ",achieved)
  end
  achieved
end

# Test if at least numsubgoals components of testgoal match a component of goal (in any order)
#=
function my_test_goal( testgoal::Goal, goal::Goal, numsubgoals::Integer )
  filter( x -> (Base.findfirst( y -> y == x, testgoal ) != nothing ), goal) 
end
=#

