# Tests insertion and deletion of gates.
using Test

function test_indel( n::Int64, p::Parameters )
  @testset "testing insertion and deletion of gates" begin
    for i = 1:n
      c = random_chromosome(p)
      sav_c = deepcopy(c)
      output_values(sav_c)  # needed to set active gates to active
      insert_gate!(c)
      @test output_values(c) == output_values(sav_c) 
      delete_gate!(c)
      @test output_values(c) == output_values(sav_c) 
    end
  end
end

function test_delete_gate( n::Int64, p::Parameters )
  for i = 1:n
    c = random_chromosome(p)
    sav_c = deepcopy(c)
    output_values(sav_c)  # needed to set active gates to active
    #dg = rand(1:(c.params.numinteriors-p.numoutputs))
    #println("dg: ",dg,"  active: ",sav_c.interiors[dg].active)
    #delete_gate!(c,dg)
    #if sav_c.interiors[dg].active == false
    delete_gate!(c)
    if number_active_gates(sav_c) < sav_c.params.numinteriors
      println("test for i = ",i)
      if output_values(c) != output_values(sav_c)
        println((output_values(c),output_values(sav_c)))
        print_circuit(c)
        print_circuit(sav_c)
        return (sav_c,c)
      end
    end
  end
end 

function test_insert_gate( n::Int64, p::Parameters )
  for i = 1:n
    c = random_chromosome(p)
    sav_c = deepcopy(c)
    insert_gate!(c)
    if output_values(c) != output_values(sav_c)
      println((output_values(c),output_values(sav_c)))
      print_circuit(c)
      print_circuit(sav_c)
      return (sav_c,c)
    end
  end
end 
