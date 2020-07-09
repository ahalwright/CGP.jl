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
  println("b4: ",b4)
  println("b3: ",b3)
  println("b2: ",b2)
  println("b1: ",b1)
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
  println("factor1: ",map(x->convert(Int64,x),factor1))
  f2 = [reverse(collect(UInt64(0x00):UInt64(0x07))) for i = 0:7]
  factor2 = vcat(f2...)
  println("factor2: ",map(x->convert(Int64,x),factor2))
  prod = map(*,factor1,factor2)
  #println("prod: ",prod)
  bb = [ map(x->convert(Int64,((0x2^i) & x)>>i),prod) for i = 5:-1:0 ]
  for i = 1:6
    println("bb[",i,"]: ",bb[i])
  end
  map(convert_binary_list_to_UInt,bb)
end

function evolve_two_bit_multiplier(iterations::Int64; num_mutations::Int64=1)
  println("num_mutations: ",num_mutations)
  max_steps = 100000
  steps_list = zeros(Int64,iterations)
  p2b = Parameters(numinputs=4,numoutputs=4,numlevelsback=16,numinteriors=10)
  mfuncs = [Func(&, 2, "&"),Func(xor, 2, "xor") ]
  #mfuncs = [Func(&, 2, "&"),Func(|, 2, "or"), Func(xor, 2, "xor") ]
  goal = goal_for_two_bit_multiplier()
  #goal = UInt16[0x40c0, 0x8000, 0x6ac0, 0xa0a0]
  #goal = UInt16[0x4c00, 0x8000, 0x6ac0, 0xa0a0]
  goallist = [goal]
  for i = 1:iterations
    rc = random_chromosome(p2b,mfuncs)
    result = mut_evolve(rc,goallist,mfuncs,max_steps,num_mutations=num_mutations)
    num_active = number_active(result[1])
    steps_list[i] = result[2]
    println("i: ",i,"  steps: ",steps_list[i],"  num active: ",num_active)
  end
  println("num_mutations: ",num_mutations)
  Statistics.mean(steps_list)
end

# Note:  MyInt must be UInt64
function evolve_three_bit_multiplier(iterations::Int64; num_mutations::Int64=1)
  println("MyInt: ",MyInt)
  println("num_mutations: ",num_mutations)
  max_steps = 50000000
  steps_list = zeros(Int64,iterations)
  p3b = Parameters(numinputs=6,numoutputs=6,numlevelsback=25,numinteriors=30)
  mfuncs = [Func(&, 2, "&"),Func(xor, 2, "xor") ]
  #mfuncs = [Func(&, 2, "&"),Func(|, 2, "or"), Func(xor, 2, "xor") ]
  goal = goal_for_three_bit_multiplier()
  #goal = UInt16[0x40c0, 0x8000, 0x6ac0, 0xa0a0]
  #goal = UInt16[0x4c00, 0x8000, 0x6ac0, 0xa0a0]
  goallist = [goal]
  for i = 1:iterations
    rc = random_chromosome(p3b,mfuncs)
    result = mut_evolve(rc,goallist,mfuncs,max_steps)
    num_active = number_active(result[1])
    steps_list[i] = result[2]
    println("i: ",i,"  steps: ",steps_list[i],"  num active: ",num_active)
  end
  println("num_mutations: ",num_mutations)
  Statistics.mean(steps_list)
end
