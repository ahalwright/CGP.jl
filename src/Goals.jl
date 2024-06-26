# Functions to implement goal functions and testing for goal functions
using Printf
export Goal, GoalList, randgoal, randgoallist, rand_env_goallist, env_goal, goal_count, goal_check, my_test_goal
export randgoal_filtered, randgoallist_filtered

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
  if numinputs <= 2
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

function randgoal(p::Parameters)
  randgoal(p.numinputs,p.numoutputs)
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

function randgoallist(ngoals::Int64, p::Parameters; repetitions::Int64=1)
  randgoallist( ngoals, p.numinputs, p.numoutputs, repetitions=repetitions )
end

# Returns a goal list of ngoals which is obtained by generating div(ngoals,repetitions)
#   goals, and then perturing these base goals repetitions-1 times to generate
#   the additional goals.
# If use_env_goals==true, then repeat goals are obtained by flipping specific bits
#    of the base goals
# If use_env_goals==false, then repeat goals are obtained by flipping random bits
#   the additional goals.
function rand_env_goallist( ngoals::Int64, numinputs::Int64, numoutputs::Int64, 
      repetitions::Int64, num_flipped_bits::Int64; perturb_goal::Bool=true )
  #println("env goallist: perturb_goal: ",perturb_goal)
  result = Vector{MyInt}[]
  if numoutputs % repetitions != 0
    error("numoutputs=",numoutputs," must be a multiple of repetitions=",repetitions," in rand_env_goallist")
  end
  #println("f_bit_list: ",f_bits_list)
  for i = 1:ngoals
    #g = rand( component_type(numinputs), div( numoutputs, repetitions ) )
    g = randgoal( numinputs, div(numoutputs,repetitions) )
    if perturb_goal
      goal = perturb_bits_goal( g, num_flipped_bits, repetitions, numinputs )   
    else
      f_bits_list = flipped_bits_list( num_flipped_bits, repetitions-1, numinputs )
      goal = env_goal( g, f_bits_list )
    end
    Base.push!(result,goal)
  end
  result
end 

# return a random goal with a minimum count of min_count in the dataframe column cdf[count_field].
# The goal is filtered by count from the dataframe column cdf[count_field].
function randgoal_filtered( cdf::DataFrame, count_field::Symbol, min_count::Int64, numinputs::Int64, numoutputs::Int64 )
  max_iters = 20    # prevent an infinite loop
  goal = randgoal( numinputs, numoutputs )
  iter = 0
  count = cdf[cdf.goal.==@sprintf("0x%x",goal[1]),count_field][1]
  while iter < max_iters && count < min_count
    goal = randgoal( numinputs, numoutputs )
    println("goal: ",goal)
    count_array = cdf[cdf.goal.==@sprintf("0x%x",goal[1]),count_field]
    println("count_array: ",count_array)
    count = count_array[1]
    iter += 1
  end
  if iter < max_iters
    (goal,count)
  else
    error("too many iterations in randgoal filtered")
  end
end

# generate a random goallist of length ngoals
# The goals are filtered by count from the dataframe cdf field count_field.
function randgoallist_filtered( cdf::DataFrame, count_field::Symbol, min_count::Int64, ngoals::Int64, numinputs::Int64, 
     numoutputs::Int64; repetitions::Int64=1)
  result = Tuple{Vector{MyInt},Int64}[]
  if ngoals % repetitions != 0
    error("ngoals=",ngoals," must be a multiple of repetitions=",repetitions," in randgoallist")
  end
  #println("  collect(1:repetitions.ngoals): ",collect(1:repetitions:ngoals))
  for i = 1:repetitions:ngoals
    (goal,count) = randgoal_filtered( cdf, count_field, min_count, numinputs, numoutputs )
    for j = 1:repetitions
      #println("(i,j): ",(i,j)," i+j-1: ",i+j-1)
      #result[i+j-1] = goal
      push!(result,(goal,count))
    end
  end
  result
end     

# The bits that will be flipped to produce the additional goals in function env_goal()
# Example (where MyInt==UInt16):
# julia> flipped_bits_list( 2,3,3)
# 3-element Array{UInt16,1}:
#  0x00c0
#  0x0030
#  0x000c
function flipped_bits_list( num_flipped_bits::Int64, repetitions::Int64, numinputs::Int64 )
  one = MyInt(1)
  if num_flipped_bits > 0
    ones = one
    for i = 1:(num_flipped_bits-1)
      ones = ones<<1 | one
    end
  else
    ones = MyInt(0)   # No flips if num_flipped_bit == 0
  end 
  fb = ones << (2^numinputs-num_flipped_bits)
  #Printf.@printf("flipped_bits_list: ones:  0x%04x  fb: 0x%04x\n",ones,fb)
  [ fb >> (num_flipped_bits*i) for i = 0:(repetitions-1) ]
end

# returns a new goal with each component of goal replicated length(flipped_bits_list) times,
#   where each replication has nbits_to_perturb bits of the original component perturbed
function env_goal( goal::Goal, flipped_bits_list::Vector{MyInt} )
  #println("env goal")
  repetitions = length(flipped_bits_list)
  new_goal = Goal()
  for c in goal
    Base.push!( new_goal, c )
    for i = 1:repetitions
      Base.push!( new_goal, xor( c, flipped_bits_list[i] ))
    end
  end
  #println("new_goal: ",new_goal)
  new_goal
end

# returns an new goal with each component of goal replicated repetitions times,
#   where each replication has nbits_to_perturb randomly chosen bits 
#   of the original component perturbed
function perturb_bits_goal( goal::Goal, nbits_to_perturb::Int64, repetitions::Int64, numinputs::Int64 )
  #println("perturb goal")
  new_goal = Goal()
  for c in goal
    push!( new_goal, c )
    for i = 1:repetitions-1
      rp = rand_perturbation( 2^numinputs, nbits_to_perturb ) 
      #Printf.@printf("rp: 0x%04x\n",rp)
      #push!( new_goal, xor( c, rand_perturbation( 2^numinputs, nbits_to_perturb )))
      push!( new_goal, xor( c, rp))
   end
  end
  #println("ptb: Goal: ",goal)
  new_goal
end

# Returns a MyInt with nnbits_to_perturb 1 bits.
# When perturbing a goal, maxlen should be 2^numinputs
function rand_perturbation( maxlength::Int64, nbits_to_perturb::Int64 )
  bits_to_perturb = rand_set( maxlength, nbits_to_perturb )
  #println("bits_to_perturb: ",bits_to_perturb)
  result = MyInt(0)
  for i in bits_to_perturb
    result |= MyInt(1) << i
  end
  result
end

# returns a random set of n_elements intgers from the integer range [0,maxlen-1]
function rand_set( maxlen::Int64, n_elements::Int64 )
  @assert n_elements < maxlen
  S = Int64[]
  while length(unique!(S)) < n_elements
    push!(S,rand(0:(maxlen-1)))
  end
  S
end

# Returns true if testgoal is equal to any of the goals in goallist defined above
function goal_check( testgoal::Goal, goallist::GoalList )
  length( filter( x -> testgoal==x, goallist )) > 0
end

# Returns the number of components of testGoal that are equal to the corresponding component of goal.
# Never used
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

# Returns a 2-tuple of:
#    1.  the  maximum number of matching subgoals of testgoal to components of goallist
#    2.  the  index of the goal of goallist that gave this maximum
# If testgoal equals all subgoals of a goal of goallist, then returns numouputs and the index of the goal in goallist
# If no component of testgoal matches the corresponding component of any goal in goallist, return 0
function goal_count( testgoal::Goal, goallist::GoalList )
   findmax( map( x -> my_test_goal( x, testgoal ), goallist ))[1] 
end

# Test if at least numsubgoals components of testgoal match a component of goal (in any order)
#=
function my_test_goal( testgoal::Goal, goal::Goal, numsubgoals::Integer )
  filter( x -> (Base.findfirst( y -> y == x, testgoal ) != nothing ), goal) 
end
=#
