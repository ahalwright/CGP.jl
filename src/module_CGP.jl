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

struct Parameters
    mu::Integer
    lambda::Integer
    mutrate::Real
    targetfitness::Real
    numinputs::Integer
    numoutputs::Integer
    nodearity::Integer
    numinteriors::Integer   # Number of gates
    numlevelsback::Integer  # For LinCircuits, number of registers
end
mutable struct LinCircuit
    # Here, numinteriors represents number of gates, numlevelsback represents number of computational registers
    # the First numoutputs registers are the output registers
    circuit_vects::Vector{Vector{MyInt}}
    params::Parameters
end
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
include("evolvable_evolvability.jl")
#include("Build_chromosome.jl")
#include("epistasis.jl")
include("Complexity.jl")
include("Degeneracy.jl")
include("Evo_dict.jl")
include("Composition.jl")
include("Robustness.jl")
include("Phenotype.jl")
include("Utilities.jl")
include("random_walk.jl")
include("Evo_dict.jl")
include("Fnc_mt.jl")
include("random_walk.jl")
include("evolvable_evolvability.jl")
end
LinCircuit = CGP.LinCircuit
using .CGP
