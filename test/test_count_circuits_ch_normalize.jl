# Test the function count_circuits_ch_normalize() which is in Chromosome.jl.
using Test
#using Main.CGP
@testset "Testing function test_count_circuits_ch_normalize()" begin

# A normalized chromosome has the inputs to every gate given in increasing order.  See normalize_chromosome() in Chromosome.jl.

# The testing technique is to generate many examples of chromosomes with these Parameters and funcs, save their circuit_ints in clist.
# Then unique(clist) is a lower bound for the result of count_circuits_ch_normalize() which should be equal if nreps is chosen to be large enough.

function test_count_circuits_ch_normalize( P::Parameters, funcs::Vector{Func}, nreps::Int64 )
  println("count_circuits_ch_normalize: ",count_circuits_ch_normalize(P,funcs))
  clist = Int128[]
  for i = 1:nreps
    rch = normalize_chromosome(random_chromosome( P, funcs ))
    cint = chromosome_to_int( rch, funcs )
    #print("cint: ",cint,"  "); print_circuit(rch)
    push!( clist, cint )
  end
  length(unique(clist))
end

P = Parameters(3,1,3,3); funcs = [default_funcs(P)[3]]
nreps = 10000
print_parameters( P )
println("nreps: ",nreps)
tccn = test_count_circuits_ch_normalize( P, funcs, nreps )
cccn = Int(count_circuits_ch_normalize(P,funcs))

@test  tccn == cccn

end # @testset
