# Tests the functions in random_walk.jl
# All tests run on 10/21/21
# These functions implement two alternate techniques for computing robustness and degree evolvability.
# The default is to store results in a goal_pair_dict
# The alternative is to store results in a goal_edge_matrix
# The test is to make sure that these give the same results when the random number generated has the same seed.
# However, this technique won't work when pmap() is used.
# So to do the test, you must change pmap to map in two places in function run_random_walks_parallel()
# Assuming that this has been done:
include("../src/random_walk.jl")
using Random
pc = Parameters(2,1,5,3)  # Chromosome parameters
pl = Parameters(2,1,5,2)  # LinCircuit Parameters
ngoals = 1
#gl = randgoal( p.numinputs, p.numoutputs )
gl = [ MyInt(i) for i = 0:2^2^p.numinputs-1 ]   # all 2-input goals
dddf=run_random_walks_parallel(4,100,gl,pc,1000,output_dict=true,save_complex=false,use_lincircuit=false)
dcdf=run_random_walks_parallel(4,100,gl,pc,1000,output_dict=true,save_complex=true,use_lincircuit=false)
mdf=run_random_walks_parallel(4,100,gl,pc,1000,output_dict=false,save_complex=false,use_lincircuit=false)
dddf=run_random_walks_parallel(4,100,gl,pl,1000,output_dict=true,save_complex=false,use_lincircuit=true)
dcdf=run_random_walks_parallel(4,100,gl,pl,1000,output_dict=true,save_complex=true,use_lincircuit=true)
mdf=run_random_walks_parallel(4,100,gl,pl,1000,output_dict=false,save_complex=false,use_lincircuit=true)

# Additional tests
nwalks = 10
steps = 1000
dd = run_random_walks( nwalks, pc, steps; output_dict=true )
M = run_random_walks( nwalks, pc, steps; output_dict=false )
ssize = 2^(2^p.numinputs)
@assert matrix_to_dict(dict_to_matrix(dd,p)) == dd
@assert dict_to_matrix(matrix_to_dict(M),p) == M

