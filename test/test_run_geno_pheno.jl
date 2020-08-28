# Simple test of function run_geno_pheno_evolution()
# include("../src/CGP.jl") # Assumes that CGP.jl is loaded
# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_run_geno_pheno_evolve.csv"
else
  csvfile = "../test/data/test_run_geno_pheno_evolve.csv"
end
allgoals=false
iterations = allgoals ? 2 : 20
numinputs = 3
numoutputs = 1
nodearity = 2
numinteriors = 8
ngoals = 1
levelsback=numinteriors-4
maxsteps = 100000
df = run_geno_pheno_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, csvfile, allgoals=allgoals )
