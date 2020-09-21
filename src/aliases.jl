export Population, Chromosome, Parameters, indiv_result_type, pop_result_type, elt_type

using Distributions
using DataFrames
using Random
using QuadGK
using CSV
using Statistics
using Distributed
using Printf
import Random.seed!
export MyInt, DIST_TYPE, Population, PopVect, IPopulation, IntRange

if !@isdefined(MyInt)  # MyInt should be defined in CGP.jl
  const MyInt = UInt8
end
if !@isdefined(MyFunc)  # MyFunc should be defined in CGP.jl
  const MyFunc = UInt128
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
const IntRange = Union{Integer, Bool, AbstractRange{Int64}, AbstractRange{Bool}}
const Goal =  Vector{MyInt}
const GoalList = Vector{Goal}
const Ones = 0xffff
#const Ones = 0xffffffff

# indiv_result_type is used to parallelize population evolution in Indiv_evolution.jl
mutable struct indiv_result_type
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  ngoals::Int64
  hamming_sel::Bool
  #avgfitness::Bool
  #robust_sel::Bool
  active_only::Bool
  maxsteps::Int64
  gl_reps::Int64
  fault_tol::Bool
  fit_limit::Float64     
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

# indiv_result_type is used to parallelize population evolution in Indiv_evolution.jl
mutable struct geno_pheno_result_type
  goallist::Vector{Vector{MyInt}}
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  ngoals::Int64
  hamming_sel::Bool
  active_only::Bool
  maxsteps::Int64
  gl_reps::Int64
  steps::Int64
  logsteps::Float64
  avgfit::Float64
  nactive::Int64
  redund::Float64
  complex::Float64
  gb_complex::Float64
  degen::Float64
  #gb_degen::Float64
  sdegen::Float64   
  f_mutinf::Float64
  mutrobust::Float64
  evolvability::Float64
end

# env_result_type is used to parallelize population evolution in Env_evolution.jl
mutable struct env_result_type
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  ngoals::Int64
  hamming_sel::Bool
  #robust_sel::Bool
  active_only::Bool
  maxsteps::Int64
  gl_reps::Int64
  num_flip_bits::Int64
  perturb_goal::Bool
  avgfitness::Bool
  perm_heuristic::Bool
  fit_limit::Float64
  steps::Int64
  avgfit::Float64   # average fitness of final chromosome
  same::Int64
  worse::Int64
  better::Int64
  nactive::Int64
  redundancy::Float64
  complexity::Float64
  degeneracy::Float64
  sdegeneracy::Float64   
end

# pop_result_type is used to parallelize population evolution in Pop_evolution.jl
mutable struct pop_result_type
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  ngoals::Int64
  popsize::Int64
  max_pop_gens::Int64
  max_indiv_steps::Int64
  hamming_sel::Bool
  robust_sel::Bool
  all_max_sel::Bool
  active_only::Bool
  maxfit::Float64
  maxrobust::Float64
  nactive::Int64
  redundancy::Float64
  complexity::Float64
  degeneracy::Float64
  sdegeneracy::Float64   
end

# pop_result_type is used to parallelize population evolution in Inf_alleles.jl
mutable struct inf_alleles_result_type
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levsback::Int64
  gl::Vector{Goal}
  popsize::Int64
  max_pop_gens::Int64
  tourn_size::Int64
  func_evals::Int64
  fitness::Float64
  gen_finished::Int64
end

mutable struct evo_result_type
  goal::Goal
  nchromes::Int64
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  maxsteps::Int64
  all_count::Int64
  evolvable_count::Int64
  evo_diff_count::Int64
end

mutable struct evo_pairs_type
  source_g::Goal
  dest_g::Goal
  numinputs::Int64
  numoutputs::Int64
  numints::Int64
  levelsback::Int64
  maxsteps::Int64
  hamming_dist::Float64
  source_g_count::Int64
  dest_g_count::Int64
  steps::Int64
end
