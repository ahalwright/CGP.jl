using DataStructures

# Finds the chromosome_ints of the circuits in the component of the genotype network of the phenotype mapped to by ch.
# Example:  
#   p = Parameters(3,1,4,4); funcs=default_funcs(p)
#   ph = [0x000d]
#   ch = pheno_evolve(p,funcs,ph,5,10_000)[1]; geno_component(p,funcs,ch)
#   Returns list of chromosome_ints
function geno_component( p::Parameters, funcs::Vector{Func}, ch::Chromosome )
  ph = output_values(ch)[1]
  chi = chromosome_to_int(ch, funcs )
  mutated = Set(Int128[chi])
  q = RBTree{Int128}();
  push!(q,chi)
  while length(q) > 0
    chi = q[1]
    delete!(q,chi)
    if !(chi in mutated)
      push!( mutated, chi )
    end
    ch = int_to_chromosome( chi, p, funcs )
    mut_list = mutate_all_neutral( ch, funcs )
    for mch in mut_list
      mchi = chromosome_to_int( mch, funcs )
      if !(mchi in mutated) && !(mchi in q)
        push!(q,mchi)
      end
    end
    #println("length(q): ",length(q),"  length(mutated): ",length(mutated))
  end
  collect(mutated)
end

# For an example, see data/2_12_23/run_geno_compsA.jl
#   and see comments below
function geno_components( p::Parameters, funcs::Vector{Func}, ph::Goal, ph_csvfile::String, field::Symbol )
  df = read_dataframe( ph_csvfile )
  @assert String(field) in names(df)
  comp_lengths = Int64[]
  phints = map(x->x-1,findall(x->x==ph[1],select(df,[field])[:,1]))
  # map(x->output_values(int_to_chromosome( x, p, funcs )), phints )
  while length(phints) > 0
    comp = map(x->Int64(x), geno_component( p, funcs, int_to_chromosome(rand(phints), p, funcs ) ) )
    #println("length(comp): ",length(comp))
    push!(comp_lengths,length(comp))
    phints = setdiff( phints, comp )
  end
  comp_lengths
end

# Summarizes component counts and component sizes in the output dataframe,
# Example:
# p = Parameters(3,1,4,4); funcs=default_funcs(p); nfuncs=length(funcs)    
# field = Symbol("ph$(p.numinputs)$(p.numinteriors)$(p.numlevelsback)$(nfuncs)")  # :ph3445
# csvfile = "../data/2_12_23/ph$(p.numinputs)$(p.numinteriors)$(p.numlevelsback)$(nfuncs).csv"  # "../data/2_12_23/ph3445.csv"
# df = run_geno_components( p, funcs, ph, csvfile, field )
function run_geno_components( p::Parameters, funcs::Vector{Func}, ph::Goal, ph_csvfile::String, field::Symbol )
# field = Symbol("ph$(numinputs)$(numinteriors)$(numlevelsback)$(nfuncs)")  comp_lengths = geno_components( p, funcs, ph, ph_csvfile, field )
  unique_comp_lengths = sort( unique( comp_lengths ), rev=true )
  dd = Dict{Int64,Tuple{Int64,Int64}}( i=>(0,0) for i in unique_comp_lengths )
  for j in sort(comp_lengths,rev=true)
    dd[j] = (dd[j][1] + 1, j)
  end
  component_tuples = [ dd[j] for j in unique_comp_lengths ]
  component_counts = map( x->x[1], component_tuples )
  component_sizes = map( x->x[2], component_tuples )
  new_df = DataFrame( :component_counts=>component_counts, :component_sizes=>component_sizes )
  new_df
end
