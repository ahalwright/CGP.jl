# Functions (like mutate) that operate on phenotypes as bit strings (or Unts).
export extract_odd_even
export half_unitation, interleave

function mutate_phenotype!( pheno::Goal, p::Parameters )
  nbits = 2^p.numinputs
  i = rand(1:length(pheno))
  bitind = rand(0:(nbits-1))
  mask = MyInt(1) << bitind
  pheno[i] = xor(pheno[i],mask)
  pheno
end
 
function mutate_all_phenotype( pheno::Goal, p::Parameters )
  nbits = 2^p.numinputs
  result_list = Goal[]
  for i = 1:length(pheno)
   for  bitind = 0:(nbits-1)
      cur_pheno = deepcopy(pheno)
      mask = MyInt(1) << bitind
      cur_pheno[i] = xor(cur_pheno[i],mask)
      push!(result_list,cur_pheno)
    end
  end
  result_list
end

function Kcomp_neighbors( pheno::Goal, kdict::Dict{UInt16, Int64}, p::Parameters )
  mut_phenos = mutate_all_phenotype( pheno, p )
  map( x->kdict[x[1]], mut_phenos )
end

function Kcomp_neighbors_stats( pheno::Goal, kdict::Dict{UInt16, Int64}, p::Parameters )
  Kcomp_nbrs = Kcomp_neighbors( pheno, kdict, p )
  [ kdict[pheno[1]], mean(Kcomp_nbrs), median(Kcomp_nbrs), std(Kcomp_nbrs), quantile(Kcomp_nbrs,0.2), quantile(Kcomp_nbrs,0.8) ]
end

# If ph is interpreted as a bit string, extract the odd bits and even bits, return as a pair of MyInts
# Example:  extract_odd_even( 0x98, 3 ) returns (0x000a, 0x0004)
function extract_odd_even( ph::Unsigned, numinputs::Int64 )
  len = 2^numinputs
  odd_result = MyInt(0)
  even_result = MyInt(0)
  odd_mask = MyInt(1) << (div(len,2)-1)
  even_mask = MyInt(1) << (div(len,2)-1)
  rmask = MyInt(1) << (len-1)
  for i = 1:len
    if i % 2 == 1
      odd_result = !iszero( rmask & ph ) ? (odd_mask | odd_result) : odd_result
      odd_mask >>= 1
    else
      even_result = !iszero( rmask & ph ) ? (even_mask | even_result) : even_result
      even_mask >>= 1
    end
    #println("i: ",i,"  rmask: ",to_binary(rmask,8),"  odd_result: ",to_binary(odd_result,8),"  even_result: ",to_binary(even_result,8))
    rmask >>= 1
  end
  (odd_result, even_result)
end

function extract_odd_even( numinputs::Int64, ph::Unsigned )
  extract_odd_even( ph, numinputs )
end

function extract_odd_even( P::Parameters, ph::MyInt )
  extract_odd_even( P.numinputs, ph )
end

# Interleaves the bits of ph1 and ph2.
# numbits should be the number of bits of ph1 and ph2
# The result will have 2*numbits bits.
# Example:  interleave( 4, 0x0004, 0x000c ) returns 0x00b0
# And extract_odd_even(8,0x00b0) returns (0x0004, 0x000c)
# Test1:  ph1 = rand(0x0000:0x000f); ph2 = rand(0x0000:0x000f); phi = interleave(4,ph1,ph2); pp=extract_odd_even(8,phi); (ph1,ph2,phi,pp)
# Test2:  pp = rand(0x0000:0x00ff);(ph1,ph2)=extract_odd_even(8,pp);phi=interleave(4,ph1,ph2); (pp,ph1,ph2,phi)
# As of 9/4/23, Test1 and Test2 don't work. Look at test/test_interleave_extract.jl.
function interleave( numbits::Int64, ph1::Unsigned, ph2::Unsigned )
  MInt = typeof(ph1)
  @assert MInt == typeof(ph2)
  result = MInt(0)
  shift = 0
  #println("A: ",to_binary(result,2*numbits))
  for i = 1:numbits
    result |= (ph1 & MyInt(1)) << shift
    #println("i: ",i,"  tb1:  ",to_binary(ph1,2*numbits))
    #println("i: ",i,"  tbr:  ",to_binary(result,2*numbits))
    ph1 >>= 1
    shift += 1
    result |= (ph2 & MyInt(1)) << shift
    #println("i: ",i,"  tb2:  ",to_binary(ph2,2*numbits))
    #println("i: ",i,"  tbr:  ",to_binary(result,2*numbits))
    ph2 >>= 1
    shift += 1
  end
  #println("E: ",to_binary(result,2*numbits))
  result
end
  
# The absolute deviation of count_ones(x) from the median of possible values of count_ones(x).
# Example with P3.numinputs == 3.  In this case, half is 4.   
#  half_unitation(0x000f, P3 ) == 0
#  half_unitation(0x00ff, P3 ) == 4
function half_unitation( x::MyInt, P::Parameters )
  half = div( 2^P.numinputs, 2 )
  abs( count_ones(x) - half )
end

# Replace by find_sub_bitstring --- obsolete
# Given an unsigned integer ph where numbits are significant, find the positions of the bit pair 10 in ph
#=
function transitions10( ph::Unsigned, numbits::Int64 )
  println("ph:       ",to_binary(ph,16))
  mask = typeof(ph)(3)   # two bits
  println("mask:     ",to_binary(mask,16))
  to_find = typeof(ph)(3)
  println("to_find:  ",to_binary(to_find,16))
  positions = Int64[]
  i = 0
  while i < numbits
    println("i: ",i,"  ph: ",to_binary(ph,16))
    if mask & ph == to_find
      push!( positions, i )
    end
    ph >>= 1
    i += 1
  end
  positions
end
=#

# Given an unsigned integer ph where numbits are significant, find the positions of the sub_bitstring  to_find  in ph.
#   bitstrings are stored as contiguous parts of unsigned integers.
# Consistent with Julia, positions are 1-based where 1 is the first (leftmost) position
# Examples:
#   julia> find_sub_bitstring( 0x8623, 2, 0x0002 )
#     3-element Vector{Int64}:
#      5
#      9
#     15
#   julia> find_sub_bitstring( 0x8673, 3, 0x0007 )
#     1-element Vector{Int64}:
#      5
function find_sub_bitstring( ph::Unsigned, mask_bits::Int64, to_find::Unsigned; numbits::Int64=16 )
  my_type = typeof(ph)
  @assert typeof(to_find) == my_type
  #println("ph:         ",to_binary(ph,numbits))
  mask = typeof(ph)(2^mask_bits-1) 
  #println("mask:       ",to_binary(mask,numbits))
  to_find = typeof(ph)(to_find)
  #println("to_find:  ",to_binary(to_find,numbits))
  positions = Int64[]
  i = 1
  while i <= numbits-mask_bits+1
    #println("i: ",i,"  ph: ",to_binary(ph,numbits))
    if mask & ph == to_find
      push!( positions, i )
    end
    ph >>= 1
    i += 1
  end
  positions
end

# Number of instances of the sub_bitstrings "01" and "10" in ph
function number_transitions( ph::Unsigned; numbits::Int64=16 )
  length(find_sub_bitstring( ph, 2, typeof(ph)(1), numbits=numbits ) ) + 
      length(find_sub_bitstring( ph, 2, typeof(ph)(2), numbits=numbits ) )
end

# Creates a dataframe with one row per phenotype and columns for ntransitions and half_unitation
function ntransitions_kdict( p::Parameters, funcs::Vector{Func}; numbits=16 )
  kdict = kolmogorov_complexity_dict(p,funcs)
  df = DataFrame( :pheno=>MyInt[], :Kcomp=>Int64[], :ntransitions=>Int64[], :half_unitation=>Int64[] )
  for ph = MyInt(0):MyInt(2^2^p.numinputs-1)
    #df_row = ( ph, kdict[ph], number_transitions( ph, 2^2^p.numinputs, numbits=2^p.numinputs ) )
    df_row = ( ph, kdict[ph], number_transitions( ph, numbits=2^p.numinputs ), half_unitation(ph,p) )
    push!(df,df_row)
  end
  df
end
