# Test function component_properties() in Fnc_mt.jl
# Verify by doing multiple calls to neutral_component()
using Test
@testset "test that all neutral mutations of neutral component are in the neutral component " begin

#include("../src/Chromosome.jl")
include("../src/Neutral_component.jl")

# nc is a neutral component (a set)
# test that all neutral mutations of circuits in nc are in nc
function test_neutral_component( p::Parameters, funcs::Vector{Func}, nc::Set; use_lincircuit::Bool=false )
  nc_list = collect(nc)
  ov = use_lincircuit ? output_values( circuit_int_to_circuit( nc_list[1], p, funcs )) : output_values(int_to_chromosome(nc_list[1],p,funcs))
  println("ov: ",ov)
  for chi in nc_list
    ch = use_lincircuit ? circuit_int_to_circuit( nc_list[1], p, funcs ) : int_to_chromosome(nc_list[1],p,funcs)
    @assert output_values(ch) == ov
    mrch = filter(cch->ov==output_values(cch), mutate_all( ch, funcs, output_outputs=false, output_circuits=true ))
    #println("mrch[1]: ",mrch[1])
    for nch in mrch
      inch = use_lincircuit ? circuit_to_circuit_int( nch, funcs ) : chromosome_to_int( nch, funcs)
      #inch = Int128(inch)
      @test inch in nc
      if inch ∉ nc
        println(" circuit int ", inch, " not in neutral component")
        return
      #else
      #  println(" circuit int ", inch, " in neutral component")
      end
    end
  end 
end

p = Parameters( 2, 1, 3, 2 )
print_parameters(p)
funcs = default_funcs(p)
rlc = rand_lcircuit(p,funcs)
ov = output_values( rlc )
nc = neutral_component( rlc, funcs );  # a set of genotype ints
println("ov: ",ov,"   length neutral component: ",length(nc))
test_neutral_component( p, funcs, nc, use_lincircuit=true )

p = Parameters( 3, 1, 3, 3 )
print_parameters(p)
funcs = default_funcs(p)
rch = random_chromosome(p,funcs)
ov = output_values( rch )
nc = neutral_component( rch, funcs );  # a set of genotype ints
println("ov: ",ov,"   length neutral component: ",length(nc))
test_neutral_component( p, funcs, nc )

end # testset
