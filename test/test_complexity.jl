# Includes alternate versions which print more information
using Combinatorics

# test_degen_complexity actually tests that the different versions of complexity agree.
function rand_X( numinputs::Int64, numints::Int64 )
  MyInt_bits = MyIntBits( MyInt )
  X = rand(collect(MyInt(0):MyInt(2^numinputs-1)),numints)
end 

function rand_tc( numinputs::Int64, numints::Int64 )
  test_complexity( rand_X( numinputs, numints ), numinputs )
end 

function test_complexity( X::Vector{MyInt},  numinputs::Int64 )
  gbX = get_bits(X,numinputs)
  println(to_binary(X,2^numinputs))
  #println("X: ",X,"  gbX: ",gbX,"  cmplx: ",complexity0(X,numinputs))
  #println("===================")
  println("X: ",X,"  gbX: ",gbX,"  cmplx: ",complexity(X,numinputs))
end
  
# Version of complexity5 which prints more information
function complexity( X::Vector{MyInt}, numinputs::Int64; base=Float64=2.0 )
  n = length(X)
  gbX = get_bits( X, numinputs )
  ent_X = entropy(gbX,base=base)
  println("X: ",X,"  gbX: ",gbX," ent_X: ",ent_X)
  for k = 1:(length(X))
    #println("k: ",k,"  ")
    #for s in combinations(X,k)
      #println("  s: ",s,"  gbs: ",get_bits(s,numinputs),"  ents: ",entropy(get_bits(s,numinputs),base=base))
    #end
  end   
  ents = [map(x->entropy(x,base=base),map(x->get_bits(x,numinputs),[s for s in combinations(X,k)])) for k = 1:(length(X))]
  #println("ents: ",ents)
  ents_avg = map(x->sum(x)/length(x), ents )
  println("ents_avg: ",ents_avg)
  summand = [ents_avg[k] - k/n*ent_X for k = 1:n]
  println("summand: ",summand)
  sum(summand)
end

# Version of complexity6 which prints more information
function complexity1( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  n = length(X)
  Xinds = collect(1:length(X))
  ssum = 0.0
  for k = 1:Int(floor(n/2)) 
    subset_pairs = [(s,setdiff(Xinds,s)) for s in combinations(Xinds,k)]
    #println("k: ",k,"  subset_pairs: ",subset_pairs)
    X_pairs = map( i->( X[subset_pairs[i][1]], X[subset_pairs[i][2]] ), collect(1:length(subset_pairs)))
    #println("X_pairs: ",X_pairs)
    gbX_pairs = [ (get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    #println("gbX_pairs: ",gbX_pairs)
    mutints = [ mutual_information(get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    #println("mutints: ",mutints)
    summand = sum(mutints)/length(mutints)
    #println("k: ",k,"  sum mutual informaiton: ",summand)
    println(summand)
    if k == Int(ceil(n/2))
      summand /= 2
    end
    ssum += summand
  end
  ssum
end

# Version of complexity7 which prints more information
function complexity0( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  n = length(X)
  IX = integration( X, numinputs, base=base )
  println("X: ",X,"  IX: ",IX)
  Xinds = collect(1:length(X))
  ssum = 0.0
  for k = 1:(n-1)
    #println("k: ",k)
    for s in combinations(Xinds,k)
      compl = X[setdiff(Xinds,s)]
      IXc = integration0(compl,numinputs,base=base)
      #println("s: ",s,"  compl: ",compl,"  IXc ",IXc)
    end
    subsets = [X[setdiff(Xinds,s)] for s in combinations(Xinds,k)]
    IXk = map( x->integration(x,numinputs,base=base), subsets )
    #println("IXk: ",IXk)
    avgIXk = sum(IXk)/length(IXk)
    println("k: ",k,"  avgIXk: ",avgIXk)
    #println("k: ",k,"  avgIXk: ",avgIXk,"  summand: ",k/n*IX - avgIXk)
    ssum += k/n*IX - avgIXk
  end
  ssum
end

# Version of integration() with get_bits() that prints more information
function integration0( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
  gb_list = [ get_bits([x],numinputs) for x in X]
  ent_list =  [ entropy(get_bits([x],numinputs)) for x in X ]
  #println("X: ",X,"  entX: ",entropy(get_bits(X,numinputs)))
  #println("gb_list: ",gb_list,"  ent_list: ",ent_list,"  sum: ",sum(ent_list))
  sum(ent_list) - entropy(get_bits(X,numinputs))
end 
    
# Version of integration() without get_bits() that prints more information
function integration0( X::Vector{Int64}; base::Float64=2.0 )
  ent_list =  [ entropy([x]) for x in X ]
  println("X: ",X,"  entX: ",entropy(X))
  println("ent_list: ",ent_list,"  sum: ",sum(ent_list))
  sum(ent_list) - entropy(X)
end 

# Integration of XI without get_bits
function integration0( XI::Vector{Vector{Int64}} ; base::Float64=2.0 )
  X = vcat( XI... )   # Combines all of the lists in XI into a long list
  println("X: ",X,"  entX: ",entropy(X))
  ent_list = [entropy(x) for x in XI]
  println("ent_list: ",ent_list,"  sum: ",sum(ent_list))
  sum( entropy(x) for x in XI ) - entropy( X)
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

# Agrees with the results of get_bits
function get_bits1( v::Vector{MyInt}, numinputs::Int64 )
  numbits = 2^numinputs
  result = zeros(MyInt,numbits)
  mask_shift = length(v)-1
  for i = 1:length(v)
    shift = numbits-1
    mask = MyInt(1) << shift 
    j = numbits-1
    while j > 1
    #for j = collect(numbits:-1:1)
      rr = (mask & v[i]) >> shift 
      result[j] |= rr << mask_shift
      shift -= 1
      mask >>= 1
      j -= 1
    end
    mask_shift -= 1
  end
  result
end
      
# [([z,x,y], get_bits([z,x,y],2),[entropy(get_bits([p],2)) for p in [z,y,x]],(fract_part(integration([z,x,y],2)+entropy(get_bits([z,x,y],2))))) for y = 0x0:0x7  for x = 0x0:0x7  for z = 0x0:0x7]
