# Uses build_chromosome() to specify some example circuits given in of Principles Evol. Design I (Miller et al 2000)
# Defines functions to compute goals for 2-bit and 3-bit multiplier circuits.
# Note that the computed goal for the 2-bit multiplier circuit does not agree with the output of the two_bit_multipler() circuit of Fig. 15b.
# Then defines functions to evolve 2 and 3 bit
#include("build_circuit.jl")
using Statistics

funcs = [AND,OR,XOR]
# Function of Figure 2 of Principles Evol. Design (Miller et al 2000)
function one_bit_adder()
  c = build_chromosome( (1,2,3), ((AND,[1,2]),(XOR,[1,2]), (AND,[5,3]), (XOR,[5,3]), (OR,[4,6])), (8,7))
  context = construct_context(c.params.numinputs)
  print_chromosome(c)
  outputs = execute_chromosome(c,context)
  println("output: ",outputs)
  outputs
end
 
# Function of Fig. 15b of Principles Evol. Design (Miller et al 2000)
function two_bit_multiplier()
  c = build_chromosome( 
     (1,2,3,4), 
     ((AND,[1,3]),  # gate 5
      (AND,[2,4]),  # gate 6 
      (AND,[2,3]),  # gate 7 
      (AND,[1,4]),  # gate 8 
      (AND,[5,6]),  # gate 9 
      (XOR,[7,8]),  # gate 10
      (AND,[7,10])), # gate 11 
     (11,9,10,6))   # outputs P3,P2,P1,P0

  context = construct_context(c.params.numinputs)
  print_chromosome(c)
  outputs = execute_chromosome(c,context)
  println("output: ",outputs)
  outputs
end

# result should agree with result of two_bit_multiplier()
function check_two_bit_multiplier()
  c = construct_context(4)
  g1 = context[1] 
  g2 = context[2] 
  g3 = context[3]
  g4 = context[4]
  g5 = g1 & g3 
  g6 = g2 & g4 
  g7 = g2 & g3 
  g8 = g1 & g4
  g9 = g5 & g6
  g10 = xor(g7,g8)
  g11 = g10 & g7
  p3 = g11 
  p2 = g9
  p1 = g10
  p0 = g6 
  [p3,p2,p1,p0]
end

function goal_for_two_bit_multiplier()
  factor1 = [0x3, 0x3, 0x3, 0x3, 0x2, 0x2, 0x2, 0x2, 0x1, 0x1, 0x1, 0x1, 0x0, 0x0, 0x0, 0x0]
  factor2 = [0x3, 0x2, 0x1, 0x0, 0x3, 0x2, 0x1, 0x0, 0x3, 0x2, 0x1, 0x0, 0x3, 0x2, 0x1, 0x0]
  prod = map(*,factor1,factor2)
  b4 = map(x->convert(Int64,(0x8&x)>>3),prod)
  b3 = map(x->convert(Int64,(0x4&x)>>2),prod)
  b2 = map(x->convert(Int64,(0x2&x)>>1),prod)
  b1 = map(x->convert(Int64,(0x1&x)),prod)
  #println("b4: ",b4)
  #println("b3: ",b3)
  #println("b2: ",b2)
  #println("b1: ",b1)
  return map(convert_binary_list_to_UInt,[b4,b3,b2,b1])
end

function convert_binary_list_to_UInt( list::Vector{Int64} )
  result = MyInt(0)
  shift = 0
  for i = length(list):-1:1
    result += convert(MyInt,list[i]) << shift
    shift += 1
  end
  result
end

function UInt_to_binary_list( x::UInt64 )
  size = 64  # If the type of x is UInt64, must be 64.  If the type of x is UInt32, must be 32.
  shift = size-1
  result = fill(0,size)
  for i = 1:size
    result[i] = (x >> shift) & convert(UInt64,1)
    shift -= 1
  end
  result
end
    
function goal_for_three_bit_multiplier()
  f1=[fill(UInt64(i),8) for i=collect(0x7:-1:0x0)] 
  factor1 = vcat(f1...)
  #println("factor1: ",map(x->convert(Int64,x),factor1))
  f2 = [reverse(collect(UInt64(0x00):UInt64(0x07))) for i = 0:7]
  factor2 = vcat(f2...)
  #println("factor2: ",map(x->convert(Int64,x),factor2))
  prod = map(*,factor1,factor2)
  #println("prod: ",prod)
  bb = [ map(x->convert(Int64,((0x2^i) & x)>>i),prod) for i = 5:-1:0 ]
  for i = 1:6
    #println("bb[",i,"]: ",bb[i])
  end
  map(convert_binary_list_to_UInt,bb)
end

# Evolves a two_bit multipler circuit with the goal defined by goal_for_two_bit_multiplier() defined above
function evolve_two_bit_multiplier( numinteriors::Int64=10, max_steps::Int64=100000; 
      iterations::Int64=1, num_mutations::Int64=1, outfile::String="")
  println("num_mutations: ",num_mutations)
  steps_list = zeros(Int64,iterations)
  p2b = Parameters(numinputs=4,numoutputs=4,numlevelsback=numinteriors,numinteriors=numinteriors)
  mfuncs = [Func(&, 2, "&"),Func(xor, 2, "xor") ]
  #mfuncs = [Func(&, 2, "&"),Func(|, 2, "or"), Func(xor, 2, "xor") ]
  current_chromosome = random_chromosome(p2b,mfuncs)  # not used but establishes rc as a top-level scope variable
  goal = goal_for_two_bit_multiplier()
  #goal = UInt16[0x40c0, 0x8000, 0x6ac0, 0xa0a0]
  #goal = UInt16[0x4c00, 0x8000, 0x6ac0, 0xa0a0]
  goallist = [goal]
  for i = 1:iterations
    rc = random_chromosome(p2b,mfuncs)
    result = mut_evolve(rc,goallist,mfuncs,max_steps,num_mutations=num_mutations)
    current_chromosome = result[1]
    num_active = number_active(result[1])
    steps_list[i] = result[2]
    println("i: ",i,"  steps: ",steps_list[i],"  num active: ",num_active)
  end
  println("max_steps: ",max_steps)
  println("numinteriors=",numinteriors)
  println("num_mutations: ",num_mutations)
  if iterations==1
    steps = steps_list[1]
  else
    steps = Statistics.mean(steps_list)
  end
  fname = length(outfile)==0 ? "../data/7_10/two_bit_mult.csv" : outfile
  num_active = number_active(current_chromosome)
  println("evolve, steps, num_active, num_muts, num_ints, max_steps")
  println("2bit_mul,", steps,",",num_active,",",num_mutations,",",numinteriors,",",max_steps)
  open(fname,"a") do f
    println(f, "evolve, steps, num_active, num_muts, num_ints, max_steps")
    println(f, "2bit_mul,", steps,",",number_active(current_chromosome),",",num_mutations,",",numinteriors,",",max_steps)
  end
  if steps < max_steps
    println("output: ",sort(execute_chromosome(current_chromosome, construct_context(p2b.numinputs))))
    println("goal:   ",sort(goal))
    @assert sort(execute_chromosome(current_chromosome, construct_context(p2b.numinputs))) == sort(goal)
  end
  current_chromosome
end

# Evolves a three_bit multipler circuit with the goal defined by goal_for_three_bit_multiplier() defined above
# Usually succeeds (based on a sample of 4) with numinteriors=35 and max_steps=50000000
# Note:  MyInt must be UInt64
function evolve_three_bit_multiplier( numinteriors::Int64, max_steps::Int64; 
      iterations::Int64=1, num_mutations::Int64=1, outfile::String="")
  @assert MyInt == UInt64
  println("num_mutations: ",num_mutations)
  steps_list = zeros(Int64,iterations)
  p3b = Parameters(numinputs=6,numoutputs=6,numlevelsback=25,numinteriors=numinteriors)
  mfuncs = [Func(&, 2, "&"),Func(xor, 2, "xor") ]
  #mfuncs = [Func(&, 2, "&"),Func(|, 2, "or"), Func(xor, 2, "xor") ]
  current_chromosome = random_chromosome(p3b,mfuncs)  # not used but establishes rc as a top-level scope variable
  goal = goal_for_three_bit_multiplier()
  goallist = [goal]
  for i = 1:iterations
    rc = random_chromosome(p3b,mfuncs)
    result = mut_evolve(rc,goallist,mfuncs,max_steps)
    current_chromosome = result[1]
    num_active = number_active(result[1])
    steps_list[i] = result[2]
    println("i: ",i,"  steps: ",steps_list[i],"  num active: ",num_active)
  end
  println("max_steps: ",max_steps)
  println("numinteriors=",numinteriors)
  println("num_mutations: ",num_mutations)
  if iterations==1
    steps = steps_list[1]
  else
    steps = Statistics.mean(steps_list)
  end
  num_active = number_active(current_chromosome)
  fname = length(outfile)==0 ? "../data/7_10/three_bit_mult.csv" : outfile
  println("evolve, steps, num_active, num_muts, num_ints, max_steps")
  println("3bit_mul,", number_active,",", steps,",",num_mutations,",",numinteriors,",",max_steps)
  open(fname,"a") do f
    println(f, "evolve, steps, num_active, num_muts, num_ints, max_steps")
    println(f, "3bit_mul,", steps,",",num_active(current_chromosome),",",num_mutations,",",numinteriors,",",max_steps)
  end
  if steps < max_steps
    println("output: ",sort(execute_chromosome(current_chromosome, construct_context(p3b.numinputs))))
    println("goal:   ",sort(goal))
    @assert sort(execute_chromosome(current_chromosome, construct_context(p3b.numinputs))) == sort(goal)
  end
  current_chromosome
end
