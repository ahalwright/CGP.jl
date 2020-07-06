# Test the inf_alleles() function from src/Evolve.jl
# to run from the test folder:  julia -p 4 -L CGP.jl test_inf_alleles.jl

#include("../src/CGP.jl")
using Main.CGP
MyInt = Main.CGP.MyInt

nreps = 4
numinputs = 3
numoutputs = 2
numints = 10
levsback = 10
popsize = 5:5:20
max_pop_gens = 50000
ngoals = 8   # number goals in goallist
tourn_size = 0  
#println("pwd: ",pwd())
# The following lines allow this file to included from CGP.jl, CGP.jl/test, and CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_inf_alleles.csv"
else
  csvfile = "../test/data/test_inf_alleles.csv"
end
#println("csvfile: ",csvfile)
#inf_alleles( p, popsize, max_pop_gens, gl, tourn_size=tourn_size )
df = run_inf_alleles(nreps,numinputs,numoutputs,numints,ngoals,levsback,popsize,max_pop_gens,tourn_size,csvfile)

