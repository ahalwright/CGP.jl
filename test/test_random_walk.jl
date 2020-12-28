# Tests the functions in random_walk.jl
# These functions implement two alternate techniques for computing robustness and degree evolvability.
# The default is to store results in a goal_pair_dict
# The alternative is to store results in a goal_pair_matrix
# The test is to make sure that these give the same results when the random number generated has the same seed.
# However, this technique won't work when pmap() is used.
# So 
