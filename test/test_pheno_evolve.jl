# Test the function run_pheno_evolve() and the two versions of pheno_evolve() in Evolve.jl.

using Test
@testset "Test the function run_pheno_evolve() and the two versions of pheno_evolve() in Evolve.jl" begin
p = Parameters(3,1,8,5)
funcs = default_funcs(p)
ph = [0x0049]
max_tries = 5
max_steps = 200_000
(c,steps) = pheno_evolve( p, funcs, ph, max_tries, max_steps )
print_circuit(c)
@test output_values(c) == ph

num_circuits_per_goal = 3
c_steps_list = pheno_evolve( p, funcs, ph, num_circuits_per_goal, max_tries, max_steps );
circuits_list = map(x->x[1], c_steps_list );
map( print_circuit, circuits_list );
steps_list = map(x->x[2], c_steps_list )
map( c->(@test output_values(c)==ph), circuits_list )

phlist = [[0x0049],[0x00b6]]
df = run_pheno_evolve( p, funcs, phlist, max_tries, max_steps )
@test size(df)[2] == 10
end # testset
