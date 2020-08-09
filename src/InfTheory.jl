# Compute robustness/degeneracy, redundancy, and complexity based on the definitions of
#   the Macia and Sole (2009) paper.  Macia's definitions are based on papers by Tononi et al.
# Note that Macia uses Tononi's definition of robustness as his definition of degeneracy.
using Combinatorics
export get_probs, get_bits, degeneracy, degeneracy1 
#include("../../information_theory/src/entropy.jl")
export degeneracy, degeneracy1, complexity4, complexity5, complexity6, complexity7, redundancy, integration
export mutinf1, mutinf2, test_MyInt, MyIntBits

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
# On May 1 replaced call to get_probs with a call to get_bits
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

function complexity6( c::Chromosome; base=Float64=2.0 )
  # Macia & Sole:  Z = n = number of "interacting units".  See comments in file macia_circuit.jl.
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  complexity6( X, numinputs; base=Float64=2.0 )
end

# Tononi complexity as defined by equation 6 of Tononi et al. (1994) and eqn. 2.9 of Macia and Sole 2009.
# Based on the sum of average mutual information between subsets of X and their complements
# See Figure 2b of Tononi Edelman and Sporns
# Version of complexity6 which prints more information
function complexity6( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  n = length(X)
  Xinds = collect(1:length(X))
  ssum = 0.0
  for k = 1:Int(floor(n/2))
    subset_pairs = [(s,setdiff(Xinds,s)) for s in combinations(Xinds,k)]
    #println("k: ",k,"  subset_pairs: ",subset_pairs)
    X_pairs = map( i->( X[subset_pairs[i][1]], X[subset_pairs[i][2]] ), collect(1:length(subset_pairs)))
    #println("X_pairs: ",X_pairs)
    #gbX_pairs = [ (get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    #println("gbX_pairs: ",gbX_pairs)
    mutints = [ mutinf(get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    #println("mutints: ",mutints)
    ssum += sum(mutints)/length(mutints) 
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
