# Run in CGP.jl/src
include("../src/CGP.jl")
using Main.CGP
#using Base.Test

const numinputs = 3
const numoutputs = 4
const nodearity = 2
const numinteriors =  4
const numlevelsback = 2
if numinputs == 2
  const context = [0xC, 0xA]
elseif numinputs == 3
  const context = [0xF0, 0xCC, 0xAA]
end


funcs = default_funcs()
p = Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)

#for _ = 1:100
    c = random_chromosome(p, funcs)
    print_chromosome( c )
    # Executing on these inputs tests all possible bit combinations for inputs
    println(execute_chromosome(c, context))
#end
