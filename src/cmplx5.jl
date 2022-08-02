# Tononi complexity as defined by equation 5 of Tononi et al. (1994).
# Note:  no use of mutual information, just entropy
function cmplx5( c::Chromosome; base=Float64=2.0 )
  # Macia & Sole:  Z = n = number of "interacting }units".  See comments in file macia_circuit.jl.
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  cmplx5( X, numinputs; base=Float64=2.0 )
end
#=
function cmplx5( X::Vector{MyInt}, numinputs::Int64; base=Float64=2.0 )
  n = length(X)
  #println("X:    ",X)
  println(to_binary_matrix( X, numinputs ))
  gbX = get_bits( X, numinputs )
  println("gbX: ",gbX)  
  println(to_binary_matrix( gbX, numinputs ))
  ent_X = CGP.entropy(gbX,base=base)
  println("ent_X: ",ent_X)
  ents = [map(x->myent(x),map(y->get_bits(y,numinputs),[s for s in combinations(X,k)])) for k = 1:(length(X))]
  println("ents: ",ents)
  ents_avg = map(x->sum(x)/length(x), ents )
  println("ents_avg: ",ents_avg)
  summand = [ents_avg[k] - k/n*ent_X for k = 1:n]
  println("summand: ",summand)
  for k = 1:length(X)
    for y in [s for s in combinations(X,k)]
      gb = get_bits(y,numinputs)
      print("k: ",k,"  s: ",y,"  get_bits: ",gb)
      myent(gb)
    end
  end
  sum(summand)
end
=#

# A version of Tononi complexity for the vector of node outputs which is based on equation 5 of the 1994 paper by Tononi et al.
# This version closely follows equation 5 but turns out to be considerably less efficient than the previous get_bits() version
#   from IntTheory.jl.
function cmplx5( X::Vector{MyInt}, numinputs::Int64; base=Float64=2.0 )
  ncombs(n::Int64,k::Int64) = factorial(n,k)/factorial(n-k)
  n = length(X)
  println("X:    ",X)
  B = to_binary_matrix( X, numinputs )
  #println("B: ",B)
  entB = my_entropy( B )
  cmplx_sum = 0.0
  for k = 1:n
    #println("k: ",k,"  combs: ",[s for s in combinations((1:size(B)[1]),k)])
    ents_avg= reduce(+, map( s->my_entropy(B,s), combinations((1:size(B)[1]),k) ))/ncombs(size(B)[1],k)
    #println("k: ",k,"  ents_sum: ",ents_avg)
    cmplx_sum += (ents_avg - (k/n)*entB)
  end
  cmplx_sum
end

function myent(x)
  entx=CGP.entropy(x)
  println("ent x: ",x,"  ent(x): ",entx)
  entx
end

# Note:  to_binary() is defined in InfTheory.jl
function to_binary_matrix( row_list::Vector{MyInt}, numinputs::Int64 )
  vcat( map(rl->to_binary( rl, 2^numinputs )', row_list )... )
end

function to_hex( V::Vector{Int64}, len::Int64 )
  result = MyInt(0)
  for i = len:-1:1
    v = V[i]
    #println("i: ",i,"  v: ",v)
    result <<= 1
    result += MyInt(v)
    #println("i: ",i,"  result: ",result)
  end
  result
end       

function binary_matrix_to_hex_columns( B::Matrix, rows::Vector{Int64}=collect(1:size(B)[1])  )
  sz = length(rows)
  myvect = map(i->to_hex(B[rows,i],sz), 1:size(B)[2] )
end

function my_entropy( B::Matrix, rows::Vector{Int64}=collect(1:size(B)[1]) )
  myvect = map(i->to_hex(B[rows,i],length(rows)),1:size(B)[2])
  CGP.entropy(myvect)
end
function my_entropy( B::Matrix, rows::Vector{MyInt}=collect(MyInt(1):MyInt(size(B)[1])) )
  myvect = map(i->to_hex(B[rows,i],length(rows)),1:size(B)[2])
  CGP.entropy(myvect)
end

function mutinf_complement( B::Matrix, rows::Vector{Int64}=collect(1:size(B)[1]) )
  my_entropy(B,rows) + my_entropy(B,setdiff(collect(1:size(B)[1]),rows)) - my_entropy(B)
end

function cmplx6( c::Chromosome; mutinf::Function=mutinf1, base=Float64=2.0 )
  # Macia & Sole:  Z = n = number of "interacting units".  See comments in file macia_circuit.jl.
  #n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  cmplx6( X, numinputs, base=Float64=2.0 )
end

# A version of Tononi complexity for the vector of node outputs which is based on equation 6 of the 1994 paper by Tononi et al.
# This version closely follows equation 4 but turns out to be considerably less efficient than the previous get_bits() version
#   from IntTheory.jl.
function cmplx6( X::Vector{MyInt}, numinputs::Int64; base=Float64=2.0 )
  n = length(X)
  println("X: ",X)
  Xinds = collect(1:n)
  B = to_binary_matrix( X, numinputs )
  #println("B: ",B)
  ssum = 0.0
  for k = 1:Int(floor(n/2))
    subsets = [s for s in combinations(Xinds,k)]
    #println("subsets: ",subsets)
    mutints = map( s->mutinf_complement( B, s ), subsets )
    #println("k: ",k,"  mutints: ",mutints)
    summand = sum(mutints)/length(mutints)
    println("k: ",k,"  sum mutual information: ",summand)
    if k == Int(ceil(n/2))
      summand /= 2
    end
    ssum += summand
  end
  ssum
end     

function test_alt_get_bits( c::Chromosome )
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  if number_active( c ) == 0   # Make sure chromosome has been executed
    execute_chromosome( c, construct_context(numinputs ) )
  end
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0
  gbX = to_binary_matrix(get_bits( X, numinputs ), numinputs )
  X = to_binary_matrix( X, numinputs )
  if !(X'==gbX[end:-1:1,(end-c.params.numinteriors+1):end])
    error( (X,gbX) )
  end
end

# Runs horizontal_chromsome() with optional output of datafram to a csvfile
function run_horizontal_chromosome( numinputsrange::UnitRange, nreps::Int64; funcs::Vector{Func}=default_funcs(numinputs[1]), csvfile::String="" )
  df = DataFrame( :numinputs=>Int64[], :cmplx6_avg=>Float64[] )
  for numinputs in numinputsrange
    cmplx6_sum = 0.0
    for i = 1:nreps
      hc = horizontal_chromosome( numinputs, funcs=funcs )
      cmplx6_sum += cmplx6( hc )
    end
    println("numinputs: ",numinputs,"  cmplx6_sum: ",cmplx6_sum)
    push!( df, (numinputs,cmplx6_sum/nreps) )
  end 
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# Results from horizontal_chromosome with 1 output per interior node")
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ",funcs)
      println(f,"# MyInt: ",MyInt)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
    
# constructs a multi-output chromosome where every interior node corresponds to an output
function horizontal_chromosome( numinputs::Int64; funcs::Vector{Func}=default_funcs(numinputs) )
  p = Parameters( numinputs, numinputs, numinputs, numinputs )
  nfuncs = length(funcs)
  inputs = [ InputNode( i, false, MyInt(0)) for i = 1:p.numinputs ]
  interiors = [ InteriorNode( funcs[rand(1:nfuncs)], [ rand(1:p.numinputs) for j = 1:p.nodearity ], false, MyInt(0) ) for i = 1:p.numinteriors ]
  outputs = [ OutputNode( i + p.numinputs ) for i = 1:p.numinteriors ]
  Chromosome( p, inputs, interiors, outputs, 0.0, 0.0 )
end  
#=
# Tononi complexity as defined by equation 6 of Tononi et al. (1994) and eqn. 2.9 of Macia and Sole 2009.
# Based on the sum of average mutual information between subsets of X and their complements
# See Figure 2b of Tononi Edelman and Sporns
function cmplx6( X::Vector{MyInt}, numinputs::Int64; mutinf::Function=mutinf1, base::Float64=2.0 )
  n = length(X)
  Xinds = collect(1:length(X))
  ssum = 0.0
  for k = 1:Int(floor(n/2))
    subset_pairs = [(s,setdiff(Xinds,s)) for s in combinations(Xinds,k)]
    println("k: ",k,"  subset_pairs: ",subset_pairs)
    X_pairs = map( i->( X[subset_pairs[i][1]], X[subset_pairs[i][2]] ), collect(1:length(subset_pairs)))
    println("X_pairs: ",X_pairs)
    gbX_pairs = [ (get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    println("gbX_pairs: ",gbX_pairs)
    mutints = [ mutinf(get_bits(Xp[1],numinputs), get_bits(Xp[2],numinputs)) for Xp in X_pairs ]
    println("k: ",k,"  mutints: ",mutints)
    summand = sum(mutints)/length(mutints)
    println("k: ",k,"  sum mutual information: ",summand)
    if k == Int(ceil(n/2))
      summand /= 2
    end
    ssum += summand
  end
  ssum
end 
=#

# Tononi complexity as defined by equaion 4 of Tononi et al. (1994).
# See Figure 2a of Tononi Edelman and Sporns
function cmplx7( c::Chromosome; base::Float64=2.0 )   
  n = c.params.numinteriors 
  numinputs = c.params.numinputs
  outputs = output_values(c)
  (IN, X, O) = node_values( c )  # lists of cached values of nodes, note letter O vs digit 0  
  cmplx7( X, numinputs )
end

# Tononi complexity as defined by equaion 4 of Tononi et al. (1994).
# See Figure 2a of Tononi Edelman and Sporns
function cmplx7( X::Vector{MyInt}, numinputs::Int64; base::Float64=2.0 )
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
