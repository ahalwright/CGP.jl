# Run after including ../src/CGP.jl
#  or uncomment the following line
include("../src/CGP.jl")

# Note that the parameters are set in the file ../src/SetParams.jl
# The default functions are defined in ../src/Func.jl
# The type  MyInt is defined in CGP.jl

#using Base.Test
using Printf
using Main.CGP


funcs = default_funcs()
const context = construct_contexts(Main.CGP.numinputs)[Main.CGP.numinputs]

# Generate random circuits and count the logic functions that are computed.
# Similar to results of Raman and Wagner 2011 to generate their Figure 2b.
#num_iterations = 50
num_iterations = 500000
println("num_iterations: ",num_iterations)
counts = create_count_function_hash( p.numinputs, p.numoutputs )
for _ = 1:num_iterations
    c = random_chromosome(p, funcs)
    #print_chromosome( c )
    # Executing on these inputs tests all possible bit combinations for inputs
    outputs = execute_chromosome(c, context)
  #=
    print("outputs: ",outputs)
    if length(outputs) > 1
      println("   hamming: ",hamming(outputs[1],outputs[2]))
    else
      println()
    end
  =#
    concat_outputs = concatenate_outputs(p.numinputs, outputs)
    #print("concat:  ")
    #Printf.@printf("0x%2x\n",concat_outputs)
    increment_count_function_hash(concat_outputs, counts )
end
if num_iterations <= 50
  print_count_function_hash(counts)
end
print_count_function_summary(counts)
