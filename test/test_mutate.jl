using Test
# Works on 5/18/21

# Tests the results of doing all mutations on a simple chromosome.
funcs = default_funcs(2)
c = circuit((1,2), ((3,XOR,2,1), (4,AND,2,2)))
@test mutate_all(c,funcs) == map(x->[MyInt(x)],
  [0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x00, 0x05, 0x05, 0x0a, 0x0a, 0x08, 0x02, 0x08, 0x02])

# Tests when insert_gate_prob==1.0 and delete_gate_prob==1.0
p = Parameters(2,1,3,2)
c = random_chromosome(p)
new_c = mutate_chromosome!(deepcopy(c),funcs,insert_gate_prob=1.0)[1]
@test num_mutate_locations(new_c,funcs) == 12
@test new_c.params.numinteriors == length(new_c.interiors)
@test length(new_c.interiors) == p.numinteriors + 1
new_c = mutate_chromosome!(deepcopy(c),funcs,delete_gate_prob=1.0)[1]
@test new_c.params.numinteriors == length(new_c.interiors)
@test length(new_c.interiors) == p.numinteriors - 1  
# c.params.numinteriors must be greater than c.param.numinputs for the above test to work.

# Shows the results of mutating a specific chromosome
# Not deterministic
function test_mutate( )
  p = Parameters(2,1,2,2)
  funcs = default_funcs(p.numinputs)
  c = circuit((1,2), ((3,XOR,2,1), (4,AND,2,2))) 
  @test num_mutate_locations(c,funcs) == 6
  print_circuit(c);
  for i = 1:num_mutate_locations( c, funcs )
    sav_c = deepcopy(c)
    mutate_chromosome!(c,funcs,i)
    print_circuit(c)
    c = sav_c
  end
end
