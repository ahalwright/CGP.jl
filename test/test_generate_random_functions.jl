# Generate random circuits and count the logic functions that are computed.
# Similar to results of Raman and Wagner 2011 to generate their Figure 2b.
# Run after including ../src/CGP.jl or uncomment the following line
#include("../src/CGP.jl")
# The important functions used are in ../src/Execute.jl

# Note that the parameters are set below
# The default functions are defined in ../src/Func.jl
# The types  MyInt and MyFunc are defined in CGP.jl
# Note that the printed output would be long, output is truncated in the functions in Execute.jl.

#using Base.Test
using Printf
using Main.CGP

# The following lines allow this file to included from CGP.jl, CGP.jl/test, or CGP.jl/src.
cwd = pwd()
if cwd[end-5:end] == "CGP.jl"
  csvfile = "$cwd/test/data/test_inf_alleles.csv"
  outfile = "$cwd/test/data/test_generate_random_functions.txt"
else
  csvfile = "../test/data/test_inf_alleles.csv"
  outfile = "../test/data/test_generate_random_functions.txt"
end
f = open(outfile,"w")
numinputs = 2
numoutputs = 2
numinteriors = 4
numlevelsback = numinteriors
funcs = Main.CGP.default_funcs(numinputs)
context = Main.CGP.construct_context(numinputs)
p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numinteriors, numlevelsback=numlevelsback )

num_iterations = 50   # Use a value <= 50 to see full untruncated output
#num_iterations = 10000
println(f,"num_iterations: ",num_iterations)
print_parameters(f,p)
counts = create_count_function_hash( p.numinputs, p.numoutputs )
for _ = 1:num_iterations
    c = random_chromosome(p, funcs)
    #print_chromosome( c )
    outputs = execute_chromosome(c, context)
    if num_iterations <= 50  
      print(f,"outputs: ",outputs)
      if length(outputs) > 1
        println(f,"   hamming: ",hamming(outputs[1],outputs[2]))
      else
        println(f)
      end
    end
    concat_outputs = concatenate_outputs(p.numinputs, outputs)
    #Printf.@printf(f,"concat_outputs: 0x%2x\n",concat_outputs)
    increment_count_function_hash(concat_outputs, counts )
    outputs
end
print_count_function_hash(f,counts)
print_count_function_summary(f,counts)
close(f)
