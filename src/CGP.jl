#=
CGP.jl is a library for using
[Cartesian Genetic Programming](http://www.cartesiangp.co.uk/) in
Julia. It is being developed at the University of Montana in Missoula,
MT for use in simulating the evolution of technology, though there is
nothing specific to that application in the library so it is (will be)
perfectly suitable for other applications as well.
=#

using Distributed
using DataFrames
#using DataFramesMeta
using StatsBase
using Combinatorics
using Printf
using Dates
using CSV
using Statistics
#using Revise
using Random

module CGP
using Distributed
using DataFrames
using StatsBase
using Combinatorics
using Printf
using Dates
using CSV
using Statistics
using Random
#const MyInt = UInt8       # Type of bit string integers used in bit functions
const MyInt = UInt16     # Type of bit string integers used in bit functions
#const MyInt = UInt32     # Type of bit string integers used in bit functions
#const MyInt = UInt64     # Type of bit string integers used in bit functions
#const MyInt = UInt128     # Type of bit string integers used in bit functions
#const MyInt = BigInt      # Type of bit string integers used in bit functions
const MyFunc = UInt128  # Type of concatenated output representation of functions
const maxints_for_degen = 20
include("aliases.jl")
include("Contexts.jl")
include("Parameters.jl")
#include("SetParams.jl")
include("Node.jl")
include("Chromosome.jl")
include("LinChromosome.jl")
include("Goals.jl")
include("Evolve.jl")
#include("Chrome_result.jl")
include("Execute.jl")
include("Func.jl")
include("Entropy.jl")
#include("Avg_mut_robustness.jl")
include("InfTheory.jl")
include("Indiv_evolution.jl")
include("Genotype_phenotype.jl")  # computations of evolvability and robustness are not ccorrect
#include("Env_evolution.jl")
include("Propsel.jl")
#include("Robust_evolve.jl")
#include("Inf_alleles.jl")
include("Analyze.jl")
include("RecordOutputs.jl")
include("Evolvability.jl")
#include("Build_chromosome.jl")
#include("epistasis.jl")
include("Complexity.jl")
include("Degeneracy.jl")
include("Evo_dict.jl")
include("Composition.jl")
#include("random_walk.jl")
end
using Main.CGP
MyInt=Main.CGP.MyInt
