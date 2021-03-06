# Revised on 4/29/21 to fix approximate comparisons with 0.0.   Works.
#include("../src/CGP.jl")
##include("../src/InfTheory.jl")
using Main.CGP
using Test

# Compare 2 floats for approximate equality where one of the floats might be zero
function approx_compare( x::Float64, y::Float64 )
  max( abs(x), abs(y) ) < 1.0e-12 || x ≈ y
end

nreps = 4
function test_degen_complexity( p::Main.CGP.Parameters, nreps::Integer )
  funcs = default_funcs(p.numinputs)
  for i = 1:nreps  
    println("i: ",i)
    c = Main.CGP.random_chromosome(p,funcs)
    #println("c: ")
    #print_chromosome(c)
    try
      @test degeneracy(c) ≈ degeneracy1(c)
    catch
      println((degeneracy(c),degeneracy1(c)))
    end
    try
      @test approx_compare( complexity5(c), complexity6(c) )
    catch
      println((complexity5(c),complexity6(c)))
    end 
    try 
      @test approx_compare( complexity5(c), complexity4(c) )
    catch 
      println((complexity5(c),complexity4(c)))
    end
    try
      @test approx_compare( complexity5(c), complexity7(c) )
    catch
      println((complexity5(c),complexity7(c)))
    end
  end
end

@testset "testing multiple versions of degeneracy and complexity" begin
#=
numinputs = 2
p = Main.CGP.Parameters( numinputs=numinputs, numoutputs=2, numinteriors=6, numlevelsback=6 )
test_degen_complexity( p, nreps )
=#
numinputs = 3
p =  Main.CGP.Parameters( numinputs=numinputs, numoutputs=3, numinteriors=7, numlevelsback=6 )
test_degen_complexity( p, nreps )
end  # @testset
