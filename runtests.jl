# This must be run from  CGP.jl.
# To run:  julia -p 4 -L src/CGP.jl runtests.jl # choose nprocs according to your computer
# As of 2/8/21 all files run, but many have no tests.
# Note that for test_random_walk.jl one test written as an assertion fails (now commented out)

println("running test/test_Entropy.jl")
@everywhere include("test/test_Entropy.jl")

println("running test/test_degen_complexity.jl")
@everywhere include("test/test_degen_complexity.jl")

println("running test/test_inf_alleles.jl")
@everywhere include("test/test_inf_alleles.jl")

println("running test/test_mut_evolve.jl")
@everywhere include("test/test_mut_evolve.jl")

println("running test/test_pop_evolve.jl")
@everywhere include("test/test_pop_evolve.jl")

println("running test/test_geno_robust_evolvability.jl")
@everywhere include("test/test_geno_robust_evolvability.jl")  # no test included

println("running test/test_generate_random_functions.jl")
@everywhere include("test/test_generate_random_functions.jl")

println("running test/test_run_geno_pheno.jl")
@everywhere include("test/test_run_geno_pheno.jl")

println("running test/test_build_chromosome.jl")
@everywhere include("test/test_build_chromosome.jl")  # no test included

println("running test/test_mutrobust_evolvability.jl")
@everywhere include("test/test_mutrobust_evolvability.jl")    # runs functions with no tests

println("running test/test_next_chromosome.jl")
@everywhere include("test/test_next_chromosome.jl")    # runs functions with no tests

println("running test/test_random_walk.jl")
@everywhere include("test/test_random_walk.jl")  # runs but assertion (commented out) fails

println("running test/test_selection.jl")
@everywhere include("test/test_selection.jl")  # no test included
