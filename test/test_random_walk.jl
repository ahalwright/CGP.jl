# Tests the functions in random_walk.jl
# These functions implement two alternate techniques for computing robustness and degree evolvability.
# The default is to store results in a goal_pair_dict
# The alternative is to store results in a goal_edge_matrix
# The test is to make sure that these give the same results when the random number generated has the same seed.
# However, this technique won't work when pmap() is used.
# So to do the test, you must change pmap to map in two places in function run_random_walks_parallel()
# Assuming that this has been done:
include("../src/random_walk.jl")
using Random
p = Parameters(3,1,7,4)
ngoals = 1
g = randgoal( p.numinputs, p.numoutputs )
Random.seed!(1); 
ddf=run_random_walks_parallel(4,100,g,p,1000,output_dict=true) 
Random.seed!(1); 
mdf=run_random_walks_parallel(4,100,g,p,1000,output_dict=false) 
#@assert ddf==mdf

# Additional tests
nwalks = 10
steps = 1000
dd = run_random_walks( nwalks, p, steps; output_dict=true )
M = run_random_walks( nwalks, p, steps; output_dict=false )
ssize = 2^(2^p.numinputs)
@assert matrix_to_dict(dict_to_matrix(dd,p)) == dd
@assert dict_to_matrix(matrix_to_dict(M),p) == M

