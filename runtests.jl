# This must be run from  CGP.jl.
# To run:  julia -p 4 -L src/CGP.jl runtests.jl

println("running test/test_Entropy.jl")
include("test/test_Entropy.jl")

println("running test/test_degen_complexity.jl")
include("test/test_degen_complexity.jl")

println("running test/test_inf_alleles.jl")
include("test/test_inf_alleles.jl")

println("running test/test_mut_evolve.jl")
include("test/test_mut_evolve.jl")

println("running test/test_pop_evolve.jl")
include("test/test_pop_evolve.jl")

println("running test/test_robust_evolve.jl")
include("test/test_robust_evolve.jl")

println("running test/test_mut_evolution.jl")
include("test/test_run_mut_evolution.jl")

println("running test/test_generate_random_functions.jl")
include("test/test_generate_random_functions.jl")
