iterations = 4
numinputs = 2:2
numoutputs = 2:2
nodearity = 2
numinteriors = 6:6
numlevelsback = 6:6
ngoals = 4:4
goallistlength=8:8
levelsback=6:6
max_pop_gens_rng = 10:10
max_indiv_steps_rng = 20:20
popsize_rng=5:5
hamming_rng = true:true
robust_rng = true:true
run_pop_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, levelsback,
    max_pop_gens_rng, max_indiv_steps_rng, popsize_rng, hamming_rng, robust_rng, "testdata.csv" )
