using Combinatorics
# Should do this via Walsh transform, but for now, just do it

# count_ones(x) is a built-in Julia function with the same purpose
function count1bits( x::MyInt, numinputs::Int64 )
  shift = (2^numinputs-1)
  mask = MyInt(1) << shift
  ssum = 0
  while mask != MyInt(0)
    ssum += (mask & x) >> shift
    #println("mask: ",mask,"  shift: ",shift,"  ssum: ",ssum)
    shift -= 1
    mask >>= 1
  end
  ssum
end

function my_ones( output::MyInt, numinputs::Int64 )
  mask = One()
  result = [ count1bits( output & MyInt(m), numinputs )/count1bits( MyInt(m), numinputs ) for m = mask:-1:1 ]
  push!(result,0.0)
  result
end
  
function my_ones( output::MyInt )
  mask = One()
  result = [ count_ones( output & MyInt(m) )/count_ones( MyInt(m) ) for m = mask:-1:1 ]
  push!(result,0.0)
  result
end

function my_mask( position::Int64 )
  one(MyInt) << (position-1)
end
  
function my_mask( positions::Vector{Int64} )
  result = zero(MyInt)
  for p in positions
    result |= my_mask( p )
  end
  result
end
  

function my_size()
  i = 1
  while my_mask(i) != zero(MyInt)
    i+=1
  end
  i-1
end 

function my_size( numinputs::Int64 )
  m = One()
  i = 1
  while (my_mask(i) & m) != zero(MyInt)
    i+=1
  end
  i-1
end 
  
function epistasis( x::Integer, p::Vector{Int64}, numinputs::Int64 )
    #println("p: ",p,"  ", [[s for s in combinations(p,k)]  for k = 1:order ] )
    x = MyInt(x)
    summand = 0
    count = 0
    for q in combinations(p,length(p)-1)
      #println("  q: ",q,"  perm_output( x, p, numinputs ): ",perm_output( x, q, numinputs ))
      summand += perm_output( x, q, numinputs )
      count += 1
    end
    #println("perm_output(x,p,numinputs): ",perm_output(x,p,numinputs),"  summand/count: ",summand/count )
    abs(perm_output(x,p,numinputs) - summand/count)
end
  
function epistasis( x::Integer, order::Int64, numinputs::Int64 )
  x = MyInt(x)
  ssum = 0
  for p in combinations( collect(1:numinputs), order )
    summand = epistasis( x, p, numinputs ) 
    ssum += summand
    #println("p: ",p,"  epistasis(p): ", summand )
  end
  ssum
end

#[ (p, my_mask(Int64(my_mask( p ))), output & my_mask(Int64(my_mask( p )))) for p in combinations(collect(1:8),2) ]

function perm_output( output::MyInt, p::Vector{Int64}, numinputs::Int64 )
  context = construct_context(numinputs)
  (output & my_mask(get_bits(context,numinputs)[my_mask(p)]+1)) != zero(MyInt) ? 1 : -1
end

function Walsh( L::Int64 )
  denom = 1.0/2.0^L
  W = zeros(Int64,L,L)
  for i = 0:(L-1)
    for j = 0:(L-1)
      W[i+1,j+1] = iseven(count_ones( i & j )) ? 1 : -1
    end
  end
  W
end

# Not correct
function Wtrans( L::Int64, x::Integer )
  result = zeros(Int64,L)
  xb = to_binary(x,L)
  for i = 0:(L-1)
    result[i+1] += xb[i+1]*(iseven( count_ones( i & xb[i+1] )) ? 1 : -1 )
  end
  result
end

# returns a list of the 1-based indices i such that count_ones(i-1) == k
function k_bit_indices( L::Int64, k::Int64 )
  findall( x->count_ones(x)==k, 0:(L-1))
end

function k_bit_epistasis( W::DenseMatrix{Int64}, k::Int64, x::MyInt )
  L = size(W)[1]
  k_inds = k_bit_indices( L, k )
  wx = W*to_binary(x,L)
  #println("k_inds: ",k_inds,"  wx: ",wx)
  sum( abs(wx[i]) for i in k_inds ) 
end

function total_epistasis( W::DenseMatrix{Int64}, x::MyInt )
  L = size(W)[1]
  wx = W*to_binary(x,L)
  sum( abs(wx[i]) for i = 2:L )
end

function reciprocal_sign_epistasis()
end

function testw( W, x )
  L=size(W)[1]
  #println( [ k_bit_epistasis(W,k,x) for k = 1:Int(log2(L)) ])
  @assert sum( k_bit_epistasis(W,k,x) for k = 1:Int(log2(L))) == total_epistasis( W,x)
end

function testw( W )
  L = size(W)[1]
  for x = zero(MyInt):MyInt(Int(log2(L)))-MyInt(1)
    testw(W,x)
  end
end

