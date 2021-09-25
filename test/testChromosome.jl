# Tests only a small subset of the functions in Chromosome.jl

# So far, tests only arity 2 gates and with maxarity=2
using Test
@testset "chromosome_to_int() and enumerate_circuits() " begin
# set p.numlevelsback to be large since enumerate_circuits() and chromosome_to_int() assume maximum
p = Parameters(2,1,3,6)   
funcs = default_funcs(p.numinputs)
nfuncs = length(funcs)
cl=enumerate_circuits(p,funcs);
cil = map(ch->chromosome_to_int(ch,funcs),cl);
@test length(cil) == length(unique(cil))
@test length(cil) == count_circuits(p,nfuncs=nfuncs)
ri = rand(cil) 
@test ri == chromosome_to_int(int_to_chromosome(ri,p))  # test that chromosome_to_int() inverts int_to_chromosome()
end # @testst
