# Recreated 12/21/12

# So far, tests only arity 2 gates and with maxarity=2
using Test
@testset "chromosome_to_int() and enumerate_circuits() " begin
p = Parameters(3,1,4,3)   
funcs = default_funcs(p.numinputs)
nfuncs = length(funcs)
cl=enumerate_circuits(p,funcs);
il = map(ch->chromosome_to_int(ch,funcs),cl);
@test length(il) == count_circuits(p,nfuncs=nfuncs)
@test il == collect(0:length(il)-1)
ri = rand(il) 
@test ri == chromosome_to_int(int_to_chromosome(ri,p))  # test that chromosome_to_int() inverts int_to_chromosome()
end # @testst
