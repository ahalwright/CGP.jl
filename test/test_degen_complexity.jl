#include("../src/CGP.jl")
##include("../src/InfTheory.jl")
using Main.CGP

nreps = 0
function test_degen_complexity( p::Main.CGP.Parameters, nreps::Integer )
  for i = 1:nreps  
    println("i: ",i)
    c = Main.CGP.random_chromosome(p,funcs)
    #println("c: ")
    #print_chromosome(c)
    try
      @assert degeneracy(c) ≈ degeneracy1(c)
    catch
      println((degeneracy(c),degeneracy1(c)))
    end
    try
      @assert complexity5(c) ≈ complexity6(c)
    catch
      println((complexity5(c),complexity6(c)))
    end
  end
end

numinputs = 2
funcs = default_funcs(numinputs)
p = Main.CGP.Parameters( numinputs=numinputs, numoutputs=2, numinteriors=6, numlevelsback=6 )
test_degen_complexity( p, nreps )
#=
numinputs = 3
p =  Main.CGP.Parameters( numinputs=numinputs, numoutputs=3, numinteriors=7, numlevelsback=6 )
funcs = default_funcs(numinputs)
test_degen_complexity( p, nreps )
=#
