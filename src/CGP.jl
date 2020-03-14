#=
CGP.jl is a library for using
[Cartesian Genetic Programming](http://www.cartesiangp.co.uk/) in
Julia. It is being developed at the University of Montana in Missoula,
MT for use in simulating the evolution of technology, though there is
nothing specific to that application in the library so it is (will be)
perfectly suitable for other applications as well.
=#


module CGP
const MyInt = UInt16   # Type of bit string integers used in bit functions
const MyFunc = UInt128  # Type of concatenated output representation of functions
include("Contexts.jl")
include("Parameters.jl")
include("SetParams.jl")
include("Node.jl")
include("Chromosome.jl")
include("Execute.jl")
include("Func.jl")

end
