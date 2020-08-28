# Compute robustness/degeneracy, redundancy, and complexity based on the definitions of
#   the Macia and Sole (2009) paper.  Macia's definitions are based on papers by Tononi et al.
# Note that Macia uses Tononi's definition of robustness as his definition of degeneracy.

using Printf
using Test
using Statistics
using Combinatorics
export get_probs, get_bits, degeneracy, degeneracy1 
#include("../../information_theory/src/entropy.jl")
export degeneracy, degeneracy1, complexity4, complexity5, complexity6, complexity7, redundancy, integration
export mutinf1, mutinf2, test_MyInt, MyIntBits, to_binary
export fmutual_information, fmi_chrome, gb_complexity, gb_complexity_chrome

# if function degeneracy() or complexity() is called with numinteriors > MyIntBits(MyInt), the results of get_bits() will overflow
#   and be inaccurate.
# maxints_for_degen is the upper limit of numinteriors for calls to degeneracy() and complexity()
function test_MyInt(max_numinteriors::Int64)
  MyInt_bits = MyIntBits( MyInt )
  if max_numinteriors > MyInt_bits && maxints_for_degen > MyInt_bits
    println("max_numinteriors: ",max_numinteriors,"  maxints_for_degen: ",maxints_for_degen,"  MyInt_bits: ",MyInt_bits)
    error("max_numinteriors > MyInt_bits && maxints_for_degen > MyInt_bits in function test_MyInt.  Run with a larger width MyInt" )
  end
end

# return the maximum number of bits in a MyInt integer
function MyIntBits( my_int::Type )
  if my_int == UInt8
    8
  elseif my_int == UInt16
    16
  elseif my_int == UInt32
    32
  elseif my_int == UInt64
    64
  elseif my_int == UInt128
    128
  else
    error("error in MyIntBits")
  end
end          

# v is a vector of outputs from a subset of a circuit
# The length of v is the number of gates in the circuit
# Each element of v is interpreted as a bitstring of length 2^numinputs
# result[i] is the average of the number of 1s at bit position i 
#    where bit positions start at 1 for the rightmost (least significant) bit position
# Example:   let v = [0xe, 0x5, 0xa] = [1110, 0101, 1010]
# get_bits(v,2) = [1/3, 2/3, 2/3, 2/3] since 
#    bit position 1 has 1 one bit and 2 zero bits,
#    bit position 2 has 2 one bits and 1 zero bits,
#    bit position 3 has 2 one bits and 1 zero bits,
#    bit position 4 has 2 one bits and 1 zero bits
# See the example below for the get_bits() function
function get_probs( v::Vector{Main.CGP.MyInt}, numinputs::Integer )
  println("get_probs(",v,")")
  result = zeros(Float64,2^numinputs)
  len_v = length(v)
  mask = MyInt(1)
  mask_shift = 0
  for i = 1:2^numinputs
    for j = 1:len_v
      result[i] += (v[j] & mask) >> mask_shift
    end
    result[i] /= len_v
    mask <<= 1
    mask_shift += 1
  end
  result
end

# v is a vector of outputs from a subset of a circuit
# The length of v is the number of gates in the circuit
# Each element of v is interpreted as a bitstring of length 2^numinputs
# The result is a vector of bitstrings of length 2^numinputs
#   where result[i][j] is 
# Example:   let v = [0xe, 0x5, 0xa] = [1110, 0101, 1010]
# get_bits(v,2) = [0x2, 0x5, 0x6, 0x5] = [0010, 0101, 0110, 0101]
# get_bits can be interpreted as a transpose operation on the bit matrix:
#    1110
#    0101
#    1010
function get_bits( v::Vector{Main.CGP.MyInt}, numinputs::Int64 )
  result = zeros(MyInt,2^numinputs)
  reverse_v = reverse(v)
  mask = MyInt(1)
  mask_shift = 0
  for i = 1:2^numinputs
    shift = 0
    for j = 1:length(v)
      result[i] |= (reverse_v[j] & mask) << (shift - mask_shift) 
      shift += 1
    end
    mask <<= 1
    mask_shift += 1
  end
  result
end

function to_binary( x::MyInt, numbits::Int64 )
  result = Int64[]
  shift = numbits-1
  mask = MyInt(1) << (numbits-1)
  for i = 1:numbits
    push!(result, (mask & x)>>shift )
    shift -= 1
    mask >>= 1
  end
  result
end

function to_binary( X::Vector{MyInt}, numbits::Int64 )
  result = zeros(Int64,length(X),numbits)
  for j in 1:length(X)
    shift = numbits-1
    mask = MyInt(1) << (numbits-1)
    for i = 1:numbits
      #push!(result, (mask & X[j])>>shift )
      result[j,i] = (mask & X[j])>>shift
      shift -= 1
      mask >>= 1
    end
  end
  result
end

# get_bits for a single MyInt
function get_bits1( x::Main.CGP.MyInt, numinputs::Int64 )
  numbits = 2^numinputs
  result = zeros(MyInt,2^numinputs)
  shift = numbits-1 
  mask = MyInt(1) << shift
  for j = collect(numbits:-1:1)
    result[j] |= (mask & x) >> shift
    shift -= 1
    mask >>= 1
  end
  result
end

# Agrees with the result of get_bits but slightly slower
function get_bits1( v::Vector{MyInt}, numinputs::Int64 )
  numbits = 2^numinputs
  result = zeros(MyInt,numbits)
  mask_shift = length(v)-1
  one = convert(MyInt,0x00000000000000001)
  for i = 1:length(v)
    shift = numbits-1
    mask = one << shift 
    j = numbits
    while j >= 1
      result[j] |= ((mask & v[i]) >> shift) << mask_shift
      shift -= 1
      mask >>= 1
      j -= 1
    end
    mask_shift -= 1
  end
  result
end
      

mutinf1(X,Y;base=Float64=2.0 ) = mutual_information(X,Y,base=base)    # Standard definition
mutinf2(X,Y;base=Float64=2.0 ) = mutual_information([X,Y],base=base)  # Sherwin definition


# degeneracy according to equation 2.4 of Macia and Sole
# On April 13 replaced Macia Z with Tononi n as the number of interacting units
# On May 1 replaced call to get_probs with a call to get_bits
# Default mutinf() function is the standard defintion
# Alternate is  mutinf=mutinf2  for the Sherwin defintion
function degeneracy( c::Chromosome; mutinf::Function=mutinf1, base=Float64=2.0 )
  # Macia & Sole:  Z = number of "interacting units".  See comments in file macia_circuit.jl.
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  gbX = get_bits( X, numinputs )
  gbO = get_bits( O, numinputs )
  mi_X_O = mutinf(gbX,gbO,base=base)
  # apply get_bits and then mutual_information to the collection of subsets produced by combinations
  mutinf_list = [map(x->mutinf(x,gbO,base=base),map(x->get_bits(x,numinputs),[s for s in combinations(X,k)]))  
        for k = 1:(length(X))]
  #println("mutinf_list: ",mutinf_list)
  mi_avg = map(x->sum(x)/length(x), mutinf_list )
  #println("mi_avg: ",mi_avg)
  summand = [mi_avg[k] - k/n*mi_X_O for k = 1:n]
  @assert summand[n] == 0.0
  #println("summand: ",summand)
  sum(summand)
end

# degeneracy according to equation 2.5 of Macia and Sole
# In order to produce sets and their complements, this function works on indices of sets
#   rather than on the sets themselves as in the defintion of degeneracy.
# On April 13 replaced Macia Z with Tononi n as the number of interacting units
# On May 1 replaced call to get_probs with a call to get_bits
function degeneracy1( c::Chromosome; mutinf::Function=mutinf1, base=Float64=2.0 )
  # Macia & Sole:  Z = number of "interacting units".  See comments in file macia_circuit.jl.
  #Z = c.params.numinteriors - c.params.numoutputs  
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(c.params.numinputs) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  gbX = get_bits( X, numinputs )
  gbO = get_bits( O, numinputs )
  mi_X_O = mutinf(gbX,gbO,base=base)
  Xinds = collect(1:length(X))
  ssum = 0.0
  for k = 1:n
    #println("k: ",k)
    # Each element of combs is a 2-tuple (s,setdiff(Xinds,s)) where s is a subset of Xinds
    #   and setdiff(Xinds,s) is its complement
    combs = [(s,setdiff(Xinds,s)) for s in combinations(Xinds,k)]
    #println("combs: ",combs)
    # convert set indices into sets.  
    # Note that if x is a list of indices, X[x] is the elements of X whose indices are in x
    Xboth=map(i->map(x->X[x],combs[i]),collect(1:length(combs)))
    #println("Xboth: ",Xboth)
    # Apply get_bits function to the sets of Xboth
    gbXboth = map(i->map(x->get_bits(x,numinputs),Xboth[i]),collect(1:length(Xboth)))
    #println("gbXboth: ",gbXboth)
    #println("length(gbXboth): ",length(gbXboth))
    miXboth = [sum(map( x->mutinf(gbXboth[i][x],gbO,base=base),[1,2]) ) - mutinf(gbX,gbO,base=base) for i =1:length(Xboth)]
    #println("miXboth: ",miXboth)
    ssum += sum(miXboth)/length(miXboth)
  end
  #println("result: ",ssum/2)
  ssum/2.0
end

# Tonini redundancy as defined by equation 2.8 of Macia and Sole (2009)
#    and by equation 3 of Tononi et al. (1999)  (eq. 2.8 is unclear)
#function redundancy( c::Chromosome; base::Float64=2.0 )
function redundancy( c::Chromosome; mutinf::Function=mutinf1, base=Float64=2.0 )
  Z = c.params.numinteriors - c.params.numoutputs
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(c.params.numinputs) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  gbX = get_bits( X, numinputs )
  gbO = get_bits( O, numinputs )
  ent_X = entropy(gbX,base=base)
  sum([ mutinf(get_bits([s],numinputs),gbO,base=base) for s in X ]) - ent_X
end

# Tononi complexity as defined by equation 5 of Tononi et al. (1994).
# Note:  no use of mutual information, just entropy
function complexity5( c::Chromosome; base=Float64=2.0 )
  # Macia & Sole:  Z = n = number of "interacting units".  See comments in file macia_circuit.jl.
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  complexity5( X, numinputs; base=Float64=2.0 )
end

function complexity5( X::Vector{MyInt}, numinputs::Int64; base=Float64=2.0 )
  n = length(X)
  gbX = get_bits( X, numinputs )
  #println("X: ",X,"  gbX: ",gbX)  
  ent_X = entropy(gbX,base=base)
  ents = [map(x->entropy(x,base=base),map(x->get_bits(x,numinputs),[s for s in combinations(X,k)])) 
        for k = 1:(length(X))]
  #println("ents: ",ents)
  ents_avg = map(x->sum(x)/length(x), ents )
  #println("ents_avg: ",ents_avg)
  summand = [ents_avg[k] - k/n*ent_X for k = 1:n]
  #println("summand: ",summand)
  sum(summand)
end

function complexity6( c::Chromosome; mutinf::Function=mutinf1, base=Float64=2.0 )
  # Macia & Sole:  Z = n = number of "interacting units".  See comments in file macia_circuit.jl.
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  complexity6( X, numinputs, mutinf=mutinf, base=Float64=2.0 )
end

# Tononi complexity as defined by equation 6 of Tononi et al. (1994) and eqn. 2.9 of Macia and Sole 2009.
# Based on the sum of average mutual information between subsets of X and their complements
# See Figure 2b of Tononi Edelman and Sporns
function complexity6( X::Vector{MyInt}, numinputs::Int64; mutinf::Function=mutinf1, base::Float64=2.0 )
  n = length(X)
  Xinds = collect(1:length(X))
  ssum = 0.0
  for k = 1:Int(floor(n/2))
  #for k = 1:n
    subset_pairs = [(s,setdiff(Xinds,s)) for s in combinations(Xinds,k)]
    println("k: ",k,"  subset_pairs: ",subset_pairs)
    X_pairs = map( i->( X[subset_pairs[i][1]], X[subset_pairs[i][2]] ), collect(1:length(subset_pairs)))
    #println("X_pairs: ",X_pairs)
    gbX_pairs = [ (get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    println("gbX_pairs: ",gbX_pairs)
    mutints = [ mutinf(get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    println("k: ",k,"  mutints: ",mutints)
    summand = sum(mutints)/length(mutints)
    println("k: ",k,"sum mutual informaiton: ",summand)
    if k == Int(ceil(n/2))
      summand /= 2
    end
    ssum += summand
  end
  ssum
end 

# Tononi complexity as defined by equaion 4 of Tononi et al. (1994).
# See Figure 2a of Tononi Edelman and Sporns
function complexity7( c::Chromosome; base::Float64=2.0 )   
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  outputs = output_values(c)
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0  
  complexity7( X, numinputs )
end

# Tononi complexity as defined by equaion 4 of Tononi et al. (1994).
# See Figure 2a of Tononi Edelman and Sporns
function complexity7( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  n = length(X)
  IX = integration( X, numinputs, base=base )
  Xinds = collect(1:length(X)) 
  ssum = 0.0
  for k = 1:(n-1)
    subsets = [X[setdiff(Xinds,s)] for s in combinations(Xinds,k)]
    IXk = map( x->integration(x,numinputs,base=base), subsets )
    avgIXk = sum(IXk)/length(IXk)
    ssum += k/n*IX - avgIXk
  end
  ssum
end
  
function integration( XI::Vector{Vector{MyInt}}, numinputs::Int64; base::Float64=2.0 )
  X = vcat( XI... )   # Combines all of the lists in XI into a long list
  sum( entropy( get_bits(x,numinputs)) for x in XI ) - entropy( get_bits(X,numinputs))
end

# Integration assuming that the components are the individual elements of X
function integration( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  ent_list =  [ entropy(get_bits([x],numinputs)) for x in X ]
  sum(ent_list) - entropy(get_bits(X,numinputs))
end

# Tononi complexity as defined by equaion 4 of Tononi et al. (1994).
# Based on the difference between the integration of X and the average integration of subsets.
# See Figure 2a of Tononi Edelman and Sporns
function complexity4( c::Chromosome; base::Float64=2.0 )
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  @assert n == length(X)
  IX = integration(X,numinputs)
  #println("IX: ",IX)
  Is = [[ integration(s,numinputs) for s in combinations(X,k) ] for k = 1:(n-1) ]
  #println("Is: ",Is)
  I_avg = map(is->sum(is)/length(is), Is)
  #println("I_avg: ",I_avg)
  #println( [ k/n*IX - I_avg[k] for k = 1:(n-1) ])
  sum( [ k/n*IX - I_avg[k] for k = 1:(n-1) ])
end

function complexity66( c::Chromosome; mutinf::Function=mutinf1, base=Float64=2.0 )
  # Macia & Sole:  Z = n = number of "interacting units".  See comments in file macia_circuit.jl.
  n = c.params.numinteriors
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(c.params.numinputs) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  gbX = get_bits( X, numinputs )
  Xinds = collect(1:length(X))
  ssum = 0.0
  #for k = 1:Int(floor(n/2))
  for k = 1:n
    println("k: ",k)
    # Each element of combs is a 2-tuple (s,setdiff(Xinds,s)) where s is a subset of Xinds
    #   and setdiff(Xinds,s) is its complement
    combs = [(s,setdiff(Xinds,s)) for s in combinations(Xinds,k)]
    #println("combs: ",combs)
    # convert set indices into sets.
    # Note that if x is a list of indices, X[x] is the elements of X whose indices are in x
    Xboth=map(i->map(x->X[x],combs[i]),collect(1:length(combs)))
    #println("Xboth: ",Xboth)
    # Apply get_bits function to the sets of Xboth
    gbXboth = map(i->map(x->get_bits(x,numinputs),Xboth[i]),collect(1:length(Xboth)))
    println("gbXboth: ",gbXboth)
    #println("length(gbXboth): ",length(gbXboth))
    miXboth = [mutinf(gbXboth[i][1],gbXboth[i][2],base=base) for i =1:length(gbXboth)]
    println("miXboth: ",miXboth)
    ssum += sum(miXboth)/length(miXboth)
  end
  ssum/2.0
end  

# The information statistics used in the paper "FRENKEN, K. and NUVOLARI, A., 2004. 
# 'The early development of the steam engine: an evolutionary interpretation using complexity theory.'
#  Industrial and Corporate Change, 13(2), pp. 419-450."
# Frenken refers to Theil 1969 and 1972 as sources.
# A limitation is that all design attributes are binary, whereas in the paper one can have higher arity attributes.
# A design is specified by a string of attributes.  In this implmentation, attributes are binary,
#   whereas in the paper attributes may be of higher arity.

# P is a vector of the designs of a collection of examples of the technology.
# Thus, P[1] is the binary string that specifies the design of example.
# marginal_distributions returns the marginal distributions of the collection of designs.
# Each marginal distribution is specified by the frequency of a 1 in that bit position
function marginal_distributions( P::Vector{MyInt}, numinputs::Int64 )
  len = length(P)
  mask = MyInt(1) << (2^numinputs-1) 
  result = zeros(Float64,2^numinputs)
  for i = 1:2^numinputs
    result[i] = sum( (x & mask != 0) ? 1 : 0 for x in P)/len
    mask >>= 1
  end
  result
end

# Returns the product of the marginal distributions corresponding to design s
# Check: P = rand(0x0000:0x00ff,19) 
#  sum( product_marginals( x, marginal_distributions(P,3)) for x = 0x0000:0x00ff ) 
#    should be approximatly 1.0
function product_marginals( s::MyInt, marginals::Vector{Float64} )
  numinputs = Int(round(log2(length(marginals))))
  result = 1.0
  mask = MyInt(1) << (2^numinputs-1) 
  for i = 1:2^numinputs
    result *= s & mask != 0 ? marginals[i] : 1.0-marginals[i]
    mask >>= 1
  end
result
end

# Based on equation (7) of Frenken and Nuvolari (2004)
# Works but assumes that all attributes are binary
function fmutual_information( P::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  #println("P: ",P)
  len = length(P)
  md =  marginal_distributions(P,numinputs)
  #println("md: ",md)
  dst = pop_to_dist(P)
  #println("length(dst): ",length(dst),"  dst: ",dst)
  ssum = 0.0
  for s in keys(dst)
    pm = product_marginals(s,md)
    summand = dst[s] > 0.0 ? dst[s]*log(base, dst[s]/product_marginals(s,md) ) : 0.0
    #@printf("s: 0x%x ",s)
    #println("  pm: ",pm,"  summand: ",summand)
    ssum += dst[s] > 0.0 ? dst[s]*log(base, dst[s]/product_marginals(s,md) ) : 0.0
  end
  ssum
end

# Test
function test_frenken( P::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  md = marginal_distributions( P, numinputs )
  dst = pop_to_dist( P )   # This collects the frequencies of various designs: each key is a desing
  println("dst: ",dst)
  pm = [ product_marginals( s, md ) for s = MyInt(0):MyInt(2^numinputs-1) ]
  println("product_marginals: ",pm)
  summands = [ (get(dst,s,0.0) > 0.0 ? dst[s]*log(base, dst[s]/product_marginals(s,md) ) : 0.0) for s = MyInt(0):MyInt(2^numinputs-1) ]
  println("summands: ",summands)
  println("frenken mutual_information: ",sum(summands))
  @test sum(summands) == fmutual_information( P, numinputs )
end
  
# The example from the file 'Frenken entropy statistics 8_17-20.docx'
# test_frenken1( MyInt[1,1,2,3], 2 )

function test_rand_frenken( nreps::Int64, num_instances::Int64, numinputs::Int64; base::Float64=2.0 )
  for i = 1:nreps
    P = randgoal( numinputs, num_instances, )
    println("P: ",P,"  fmi: ",fmutual_information( P, numinputs, base=base ))
  end
end

# x can be generated by randgoal( numinputs, 2^numinputs )
function gb_complexity( x::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  complexity5( get_bits( x, numinputs ), numinputs, base=base )
end

# Returns fmutual_information( node values of the interior nodes of c )
function gb_complexity_chrome( c::Chromosome; base::Float64=2.0 )
  output_values(c)
  gb_complexity( node_values(c)[2], c.params.numinputs, base=base)
end

# Returns fmutual_information( node values of the interior nodes of c )
function fmi_chrome( c::Chromosome )
  output_values(c)
  fmutual_information( node_values(c)[2], c.params.numinputs )
end

# x can be generated by randgoal( numinputs, 2^numinputs )
function gb_fmi( x::Vector{MyInt}, numinputs::Int64 )
  fmutual_information( get_bits( x, numinputs ), numinputs )
end

# Test the correlation between functions applied to Vectors of MyInts
function test_correlation( nreps::Int64, function1::Function, function2::Function, 
    numinputs::Int64, numoutputs::Int64 )
  result1 = zeros(Float64,nreps)
  result2 = zeros(Float64,nreps)
  for i = 1:nreps
    rand_vec = randgoal(numinputs,numoutputs)
    result1[i] = function1( rand_vec, numinputs )
    result2[i] = function2( rand_vec, numinputs )
  end
  Statistics.cor( result1, result2 )
end
