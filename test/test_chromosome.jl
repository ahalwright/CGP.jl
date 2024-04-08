# Recreated 1/13/23, works 9/9/23

# Tests only arity 2 gates and with maxarity=2
using Test
@testset "chromosome_to_int() and enumerate_circuits() " begin
p = Parameters(3,1,3,3)  # Any legal parameters should work, except that computation of il may take too long.
@assert p.numoutputs==1  # chromosome_to_int and int_to_chromosome assume 1 output
funcs = default_funcs(p)
nfuncs = length(funcs)
cl=enumerate_circuits_ch(p,funcs);
@test length(cl) == count_circuits_ch(p,nfuncs=nfuncs)
il = map(ch->chromosome_to_int(ch,funcs),cl);
@test il == collect(0:length(il)-1)
ri = rand(il) 
@test ri == chromosome_to_int(int_to_chromosome(ri,p))  # test that chromosome_to_int() inverts int_to_chromosome()
# Unfortunately, it is difficult to test int_to_chromosome(chromosome_to_int(ch) == ch because equality for chromosomes doesn't work
end # @testset

# Check that all mutations are unique
numinputs = 2
p = Parameters(numinputs,1,3,4)
funcs = default_funcs(p.numinputs)
Random.seed!(1); 
rch = random_chromosome( p, funcs )
ma = mutate_all( rch, funcs, output_circuits=true, output_outputs=false)
# @test length(ma) == num_mutate_locations(rch,funcs)  # fails 3/6/24 with length(ma) > num_mutate_locations(rch,funcs)
ima = map( ch->chromosome_to_int(ch,funcs), ma )
@test length(unique(ima)) == length(ima)   # Check that all mutations are unique

p = Parameters(2,1,3,4)
funcs = default_funcs(p)
rch=random_chromosome(p,funcs)
println("rch: ")
print_circuit(rch)
#test_remove_inactive(rch)
@test output_values(rch) == output_values( remove_inactive( rch ) )
@test count_circuits_ch( p, funcs ) == length(enumerate_circuits_ch( p, funcs ))  # 72000

end #testset
