# Test the modularity hypothesis of Kashtan 2005.
# the modular goals of Kashtan are saved in notes/diary12_15.txt

function test_switching_goals( ngens::Int64, gens_to_switch::Int64, numinteriors::Int64, levels_back::Int64,
    max_steps::Int64, max_steps_for_switch::Int64 )
  p = Parameters( 4, 1, numinteriors, levels_back )
  goal1 = [0x0660]
  goal2 = [0x6ff6]
  
  funcs = default_funcs( p.numinputs )
  c = random_chromosome( p, funcs )
  res = mut_evolve( random_chromosome( p, funcs ), [goal1], funcs, max_steps )
  println("ov: ",output_values(res[1]))
  while output_values(res[1]) != goal1
    res = mut_evolve( random_chromosome( p, funcs ), [goal1], funcs, max_steps )
    println("ov: ",output_values(res[1]))
  end
  c = res[1]
  print_build_chromosome(c)
  goal = goal2
  for gen = 1:ngens 
    print("gen: ",gen,"  ")
    res = mut_evolve( c, [goal], funcs, max_steps_for_switch )
    c = res[1]
    #println("gen: ",gen,"  ov: ",output_values(c))
    goal = (goal == goal1) ? goal2 : goal1
  end
end

function modular_goals()
  #Figure 2a of Kashtan et al. 2005
  rc=build_chromosome((1,2,3,4),(
  (NAND,[1,2]), #5
  (NAND,[3,4]), #6
  (NAND,[1,6]), #7
  (NAND,[5,3]), #8
  (NAND,[2,6]), #9
  (NAND,[5,4]), #10
  (NAND,[8,10]), #11
  (NAND,[7,9]), #12
  (NAND,[11,12]), #13
  (NAND,[13,13])), #14
  (14,))
  
  #Figure 2b epoch 1 of Kashtan et al. 2005  
  e1c=build_chromosome((1,2,3,4),(
  (NAND,[1,2]), #5
  (NAND,[1,5]), #6
  (NAND,[5,2]), #7
  (NAND,[6,7]), #8
  (NAND,[3,4]), #9
  (NAND,[3,9]), #10
  (NAND,[9,4]), #11
  (NAND,[10,11]), #12
  (NAND,[8,12]), #13
  (NAND,[8,12]), #14
  (NAND,[13,14])), #15
  (15,))   
  
  #Figure 2b epoch 2 of Kashtan et al. 2005
  e2c=build_chromosome((1,2,3,4),(
  (NAND,[1,2]), #5
  (NAND,[1,5]), #6
  (NAND,[5,2]), #7
  (NAND,[6,7]), #8
  (NAND,[3,4]), #9
  (NAND,[3,9]), #10
  (NAND,[9,4]), #11
  (NAND,[10,11]), #12
  (NAND,[8,8]), #13
  (NAND,[12,12]), #14
  (NAND,[13,14])), #15
  (15,))         
end
