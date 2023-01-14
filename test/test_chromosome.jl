# Recreated 1/13/23

# So far, tests only arity 2 gates and with maxarity=2
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
end # @testset
