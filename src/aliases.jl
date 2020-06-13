export Population 

using Distributions
using DataFrames
using Random
using QuadGK
using Distributed
import Random.seed!

if !@isdefined(MyInt)
  const MyInt = UInt8
end
const elt_type = Union{Int64,MyInt,String,Float64}
const DIST_TYPE = Dict{Any,Float64}
#const IDIST_TYPE = Dict{Int64,Float64}
#const SDIST_TYPE = Dict{String,Float64}
const Population = Union{Vector{Int64},Vector{MyInt},Vector{String},Vector{Float64}}
const PopVect = Union{Vector{Vector{Int64}},Vector{Vector{MyInt}},Vector{Vector{String}},Vector{Vector{Float64}}} 
#const IPopulation = Array{Int64,1}
#const MIPopulation = Array{MyInt,1}
#const SPopulation = Array{String,1}
const FPopulation = Array{Float64,1}
#const CPopulation = Vector{Chromosome}   # population of chromosomes

mutable struct run_result_type
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  ngoals::Int64
  hamming_sel::Bool
  #robust_sel::Bool
  active_only::Bool
  maxsteps::Int64
  steps::Int64
  same::Int64
  worse::Int64
  better::Int64
  nactive::Int64
  redundancy::Float64
  complexity::Float64
  degeneracy::Float64
  sdegeneracy::Float64   
end

