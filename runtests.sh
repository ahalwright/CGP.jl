#!/bin/bash
# Bash scrpt to run basic tests
cd test
julia -L ../src/CGP.jl test_Entropy.jl
julia -L ../src/CGP.jl test_degen_complexity.jl
julia -L ../src/CGP.jl test_inf_alleles.jl
julia -L ../src/CGP.jl test_mut_evolve.jl
julia -L ../src/CGP.jl test_pop_evolve.jl > data/pop_evolve.txt
julia -L ../src/CGP.jl test_robust_evolve.jl
julia -L ../src/CGP.jl test_run_mut_evolution.jl > data/run_mut_evolution.txt
julia -L ../src/CGP.jl test_generate_random_functions.jl > data/test_generate_random_functions.txt
cd ..

