# Basic setup to do tests

include("../src/CGP.jl")
using Main.CGP
funcs = default_funcs()
context = construct_contexts(p.numinputs)[p.numinputs]
c = random_chromosome( p, funcs )
execute_chromosome(c,context)
print_chromosome(c)
println("num_mutate_locations: ",num_mutate_locations(c,funcs))
println("num_act: ",number_active(c))
println("num_act_old: ",number_active_old(c))
@assert number_active(c) == number_active_old(c)
mutate_chromosome!(c,funcs)
execute_chromosome(c,context)
print_chromosome(c)
println("num_mutate_locations: ",num_mutate_locations(c,funcs))
println("num_act: ",number_active(c))
println("num_act_old: ",number_active_old(c))
@assert number_active(c) == number_active_old(c)
nv = node_values( c )
println("nv: ",nv)
gb = get_bits( nv[2], p.numinputs )
println("gb: ",gb)


