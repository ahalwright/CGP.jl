# This file can be run by the following command run in the folder CGP.jl/test/
# test>  julia -p 4 -L ../src/CGP.jl runtests.jl
println("test/runtests.jl")
include("test_Entropy.jl")
include("test_degen_complexity.jl")
include("test_inf_alleles.jl")
include("test_mut_evolve.jl")
include("test_pop_evolve.jl")
include("test_robust_evolve.jl")
include("test_run_mut_evolution.jl")
include("test_generate_random_functions.jl")

