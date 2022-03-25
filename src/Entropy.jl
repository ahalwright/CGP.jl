# Extensive tests are in ../test/test_Entropy.jl
# TODO:  Write up a high-level description of the functions in this file
#using CSV
using DataFrames
export dist_check, pop_to_dist, pops_to_dist, pops_to_tbl, pop_counts_to_tbl, table_row_to_dist
export entropy, relative_entropy, rel_entropy, conditional_entropy, joint_entropy, mutual_information
export sherwin_mutual_information, row_marginal, column_marginal, adami_complexity
#include("../../information_theory/src/aliases.jl")

 
function dist_check( p::DIST_TYPE )
  sum = 0.0
  for x in keys(p)
    if p[x] < 0.0
      error(" value of distribution ", p, " is negative ")
    end
    sum += p[x]
  end
  if !isapprox(sum,1.0)
      error(" sum of distribution ", p, " is not 1.0")
  end
end

# Converts a population to a distribution (i. e., a Dict)
function pop_to_dist( pop::AbstractVector)
  N = length(pop)
  result = DIST_TYPE()
  for i in pop
    #pval = get(result,i,0.0)
    result[i] = get(result,i,0.0) + 1.0/N
    #println("i: ",i,"  pval: ",pval,"  result: ",result)
  end
  dist_check(result)
  result
end

# Assumes that both populations are indexed over the same set.
# Thus, both population must be the same length
# There is one entry in the table for each pair (P1[i],P2[i]).
# This means that the populations are ordered.
function pops_to_dist( P1::Population, P2::Population )
  len1 = length(P1)
  len2 = length(P2)
  @assert( len1 == len2 )
  result = DIST_TYPE()
  for i = 1:len1
    entry = (convert(elt_type,P1[i]),convert(elt_type,P2[i]))
    result[entry] = get(result,entry,0.0) + 1.0/len1
  end
  result
end 

# D is a DataFrame with 2 columns per population.
# Each odd numbered column is a list of Strings which are "allele names".
# The next even numbered column is an integer count for the corresponding allele name.
# Each allele name produces a column of the table.
# Example:
# julia> df = DataFrame(names1=["A","B","C"], counts1=[2,5,3], names2=["A","B","D"], counts2=[4,7,1])
# julia> pop_counts_to_tbl(df)
# 2×4 Array{Float64,2}:
#   0.0909091  0.227273  0.136364  0.0
#   0.181818   0.318182  0.0       0.0454545
function pop_counts_to_tbl( D::DataFrame, columns::Vector{Int64}=collect(1:(div(size(D)[2],2))))
  #m = div(length(names(D)),2)   # the number of rows in result
  m = length(columns)
  #println("m: ",m)
  # Convert columns of data frame D from CSV columns to columns of Strings or columns of Int64.
  for i = 1:m
    D[!,2*i-1] = convert(Vector{typeof(D[!,1][1])} ,D[!,2*i-1])
    D[!,2*i] = convert(Vector{Int64},D[!,2*i])
  end
  combined_pop = reduce(vcat, D[!,2*i-1] for i in columns)   # combine allele names
  combined_counts = reduce(vcat, D[!,2*i] for i in columns)  # combine counts
  N = length(combined_pop)
  combined_sum = sum(combined_counts)
  allele_list = unique(combined_pop)
  n = length(allele_list)   # The number columns in result
  #println("m: ",m,"  N: ",N,"  n: ",n,"  combined_sum: ",combined_sum)
  result = zeros(Float64,m,n)
  for i = 1:m
    ii = columns[i]
    #println("i: ",i,"  ii: ",ii)
    dist = Dict( allele_list .=> zeros(Float64,n))
    for j = 1:length(D[!,2*ii-1])
      dist[D[!,2*ii-1][j]] += D[!,2*ii][j]/combined_sum
    end
    result[i,:] = [ dist[a] for a in allele_list ]
  end
  result
end

# Converts a row of a joint probability distribution table to a DIST
# The DIST is normalized so that the probabilities sum to 1.0
function table_row_to_dist( tbl::Array{Float64,2}, row_index::Int64 )
  result = DIST_TYPE()
  (m,n) = size(tbl)
  ssum = sum(tbl[row_index,:])
  for j = 1:n
    result[j] = get(result,j,0.0) + tbl[row_index,j]/ssum
  end
  result
end

function row_marginal( tbl::Array{Float64,2} )
  [ sum(tbl[i,:]) for i = 1:size(tbl)[1]]
end

function column_marginal( tbl::Array{Float64,2} )
  [ sum(tbl[:,j]) for j = 1:size(tbl)[2]]
end                           

function entropy( p::DIST_TYPE; base::Float64=2.0 )
  result = 0.0
  for x in keys(p)
    result += p[x] > 0.0 ? -p[x]*log(base, p[x] ) : 0.0
    #println("x: ",x, "result_inc: ", p[x] * log(base, p[x]),"  result: ",result)
  end
  result
end 

function entropy( p::Population; base::Float64=2.0 )
  entropy( pop_to_dist(p), base=base )
end

# Entropy of a sequence of probabilities.
# entropy(get_probs( v, numinputs))  computes the entropy of a vector v
#   where each v[i] is interpeted as a bit string of length 2^numinputs
function entropy( p::Vector{Float64}; base::Float64=2.0 )
  #println("bit prob entropy")
  entf(x) = x > 0.0 ? -x*log(base,x) : 0.0
  reduce(+,map(entf,p))
end

function entropy( tbl::Array{Float64,2}, row_index::Int64; base::Float64=2.0 )
  entropy( table_row_to_dist( tbl, row_index ), base = base )
end

function entropy( x::MyInt, numinputs::Int64; base::Float64=2.0 )
  len = 2^numinputs
  cnt_ones = count_ones(x)
  entropy( [cnt_ones/len, (len-cnt_ones)/len ], base=base )
end

# relative_entropy()
# Also Kullback Leibler divergence D( q || p )
# Note that if there are any keys x of p such that q[x] == 0, the result will be NaN
#   indicating that the result is not defined.
# Reversed p and q to agree with Cover & Thomas on 9/16/19
function relative_entropy( p::DIST_TYPE, q::DIST_TYPE; base::Float64=2.0 )
  result = 0.0
  for x in keys(p)  # Assume that any keys of q that aren't in p have value 0 for p
    try qval = q[x]
    catch
      return NaN
    end
    #pval = get(q,x,0.0)
    pval = get(p,x,0.0)
    result += pval > 0.0 ? pval*(log(base,pval)-log(base,q[x])) : 0.0
    #println("x: ",x,"  pval: ",pval,"  q[x]: ",q[x],"  result: ",result)
  end
  result
end

function relative_entropy( q::Population, p::Population; base::Float64=2.0 )
  relative_entropy( pop_to_dist(q), pop_to_dist(p), base=base )
end

function relative_entropy( q::Vector{Float64}, p::Vector{Float64}; base::Float64=2.0 )
  @assert length(q) == length(p)
  result = 0.0
  for i = 1:length(p)
    if q[i] == 0.0
      return NaN
    end
    result += p[i] == 0.0 ? 0.0 : p[i]*log(base,p[i]/q[i])
  end
  result
end

# relative entropy of two incidence vectors:
function rel_entropy( q::Vector{Int64}, p::Vector{Int64} )
  q = q/sum(q)
  p = p/sum(p)
  relative_entropy( q, p )
end

# Conditional entropy of a table.
# H(X|Y)  where  X  is the column-marginal of the table, and Y is the row marginal of the table
# The rows argument specifies the rows of the table that are used.  Order of rows is not important.
# The default is all rows.
# Verification in the file /home/evotech/information_theory/test/test_cover_mutint.jl
function conditional_entropy( tbl::Array{Float64,2}, rows::Vector{Int64}=collect(1:size(tbl)[1]); base::Float64=2.0 )
  #println("rows: ",rows)
  row_sums = map(x->sum(tbl[x,:]),rows)
  tbl_sum = sum(row_sums)   # sum of the part of the table that we are using.
  result = 0.0
  for x = 1:length(rows)
    #sum_px = sum(tbl[rows[x],:])
    for y = 1:size(tbl)[2]
      pyx = tbl[rows[x],y]/row_sums[x]
      result -= tbl[rows[x],y]/tbl_sum > 0.0 ? tbl[rows[x],y]/tbl_sum*log(base,pyx) : 0.0
      #println("(x,y): ",(x,y),"  pyx: ",pyx,"  inc: ",tbl[rows[x],y]>0.0 ? -tbl[rows[x],y]*log(2,pyx) : 0.0)
    end
  end
  result
end

# joint entropy by definition of Cover and Thomas equation 2.8
function joint_entropy( tbl::Array{Float64,2}; base::Float64=2.0 )
  if ! reduce( &, tbl .>= 0.0 )
    error("all entries of tbl must be non-netative in function joint entropy")
  end
  if !(sum(tbl) ≈ 1.0)
    tbl = tbl/sum(tbl)   # Normalize to sum 1
  end
  (m,n) = size(tbl)
  entf(x) = x != 0 ? -x*log(base,x) : 0.0
  reduce(+, map(entf,tbl))
end

# Assumes that both populations are indexed over the same set.
# Thus, both population must be the same length
# There is one entry in the table for each pair (P1[i],P2[j]).
function joint_entropy( P1::Population, P2::Population; base::Float64=2.0 )
  @assert length(P1) == length(P2)
  entropy(pops_to_dist(P1,P2),base=base)
end      

# Mutual information by the defintion   I(P1;P2) = H(P1) + H(P2) - H(P1,P2) (eq. 2.45 of Cover)
# Uses the pops_to_dist() method in joint_entropy()
# Read the comments on the pops_to_dist() function:  Populations are ordered and must be the same size.
# This is the "standard" (non-Sherwin) definition of mutual information
function mutual_information( P1::Population, P2::Population; base::Float64=2.0 )
  @assert length(P1) == length(P2)
  entropy(P1,base=base) + entropy(P2,base=base) - joint_entropy(P1,P2,base=base)
end

# Mutual information by the "standard" defintion of eq. 2.28 of Cover and Thomas.
#   This is the mutual information of a joint distribution.
function mutual_information( tbl::Array{Float64,2}; base::Float64=2.0 )
  if ! reduce( &, tbl .>= 0.0 )
    error("all entries of tbl must be non-negative in function mutual information")
  end
  if !(sum(tbl) ≈ 1.0)
    tbl = tbl/sum(tbl)   # Normalize to sum 1
  end
  (m,n) = size(tbl)
  mif(i,j) = tbl[i,j] != 0.0 ? tbl[i,j]*log(base,tbl[i,j]/mr1[i]/mr2[j]) : 0.0
  mr1 = [ sum( tbl[i,:]) for i = 1:m ]
  mr2 = [ sum( tbl[:,j]) for j = 1:n ]
  reduce( +, [ mif(i,j) for i=1:m, j=1:n ] )
end

# tbl is the probability distribution for a total population which is made up of m subpopulations
#    all over the same allele set of size n.
# Row i of tbl (normalized to sum to 1) corresponds to a probability distribution over subpopulation i.
function sherwin_mutual_information( tbl::Array{Float64,2}; base::Float64=2.0 )
  @assert isapprox( sum(tbl), 1.0 )
  entf(x) = x != 0 ? -x*log(base,x) : 0.0
  (m,n) = size(tbl)
  row_weights = row_marginal(tbl)
  row_ents = [ sum( map(entf,tbl[i,:]/row_weights[i])) for i = 1:m ] # entropies of normalized rows
  #println("row_ents: ",row_ents)
  col_ents = map(entf,[sum(tbl[:,j]) for j = 1:n])   # entropies of the sums of columns
  #println("col_ents: ",col_ents)
  ent = sum(col_ents) - sum(row_ents[i]*row_weights[i] for i = 1:m)
  #println("ent: ",ent)
  ent
end

function mutual_information( D::DataFrame, columns::Vector{Int64}=collect(1:(div(size(D)[2],2))); base::Float64=2.0 )
  mutual_information(pop_counts_to_tbl( D, columns ), base=base )
end

function mutual_information( P::Population; base::Float64=2.0 )
  sherwin_mutual_information( pops_to_tbl( P ), base=base )
end

function mutual_information( P1::FPopulation, P2::FPopulation; base::Float64=2.0 )
  #println("MI pops_to_tbl")
  sherwin_mutual_information( pops_to_tbl( [P1, P2] ), base=base )
end 

function mutual_information( P::PopVect; base::Float64=2.0 )
  #println("MI pops_to_tbl")
  sherwin_mutual_information( pops_to_tbl( P ), base=base )
end 

# The Adami (physical) complexity of a population of phenotypes relative to maximum complexity
function adami_complexity( pop::Vector{MyInt}, p::Parameters )
  Hmax = p.numoutputs*2^p.numinputs  # Entropy of a population of all possbile phenotype strings
  Hmax - entropy( pop )   # Equation (5) of Adami (2002) "What is Complexity"
end
  
