# Functions to implement goal functions and testing for goal functions

export Goal, GoalList, goal_count, goal_check, my_test_goal

const Goal =  Vector{MyInt}
const GoalList = Vector{Goal}

#=
goal1 = MyInt[ 0xb, 0xe ]
goal2 = MyInt[ 0x6, 0x9 ]
goal3 = MyInt[ 0x1, 0xe ]
goal4 = MyInt[ 0x1, 0xe ]
goal5 = MyInt[ 0x3, 0x9 ]
goallist2 = [goal1, goal2, goal3 ]

goallist3 = [0x50, 0x05, 0xf5, 0xa0, 0xe0, 0x0e, 0xf1, 0x80, 0x08, 0xfa]
=#

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
function my_test_goal( testGoal::Goal, goal::Goal )
  length( filter( x -> x == 0, testGoal .⊻ goal ))
end

# Test if at least numsubgoals components of testgoal match a component of goal (in any order)
#=
function my_test_goal( testgoal::Goal, goal::Goal, numsubgoals::Integer )
  filter( x -> (Base.findfirst( y -> y == x, testgoal ) != nothing ), goal) 
end
=#

