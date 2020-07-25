# A heursitic to solve the aassignment probmem.
export match_score_perm_heuristic, perm
# 
function match_score_perm_heuristic( output::Vector{MyInt}, goal::Vector{MyInt},numinputs::Int64;
      avgfitness::Bool=false)
  #println("match_score: avgfitness: ",avgfitness)
  nc = length(output)  # number of components
  #println("output: ",output,"  goal: ",goal)
  @assert nc == length(goal)
  H = [ hamming_distance(output[i],goal[j],numinputs) for i = 1:nc, j=1:nc]
  #println("H: ",H)
  #P = permutations(collect(1:nc))
  if avgfitness
    perm( 1.0.-H )[1]
    #mxscores = [ (sum( 1.0-H[i,p[i]] for i = 1:nc ),[p[i] for i = 1:nc]) for p in P ]
    #println("mxscores: ",mxscores)
    #findmaxall(mxscores)[1][1]
  else # number of exact matches plus Hamming score for best partial match
    error("function match_score_perm_heuristic:  avgfitness must be true.")
    perm( H )
    #mxscores = [ maxscore([H[i,p[i]] for i = 1:nc]) for p in P ]
    #findmaxall(mxscores)[1]
  end
end

function rand_array( n::Int64, maxval::Int64 )
  result = rand(collect(0:maxval),n,n)/maxval
end

function perm( H::Array{Float64,2} )
  n = size( H )[1]
  assign = Dict{Int64,Int64}(i=>0 for i = 1:n) 
  notchosen = not_chosen( assign )
  for i = 1:n
    rowmaxlist = [ findmax( [ H[i,j] for j in not_chosen_cols(assign) ] ) for i in notchosen ]
    #println("rowmaxlist: ",rowmaxlist)
    amax = findmaxrand( rowmaxlist )
    #println("amax: ",amax)
    row = notchosen[amax[2]]
    ncc = not_chosen_cols(assign)
    col = ncc[amax[1][2]]
    #println("row: ",row,"  col: ",col,"  H[row,col]: ",H[row,col])
    assign[row] = col
    notchosen = not_chosen( assign )
    #println("notchosen: ",notchosen,"  not_chosen_cols: ",not_chosen_cols(assign))
  end
  p = [ assign[k] for k = 1:n  ]
  @assert allunique(p)
  fitsum = sum( H[i,p[i]] for i = 1:n )
  (fitsum,p)
end


function not_chosen( assign::Dict{Int64,Int64} )
  n = length( assign )
  filter(i->(assign[i]==0),collect(1:n))
end

function not_chosen_cols( assign::Dict{Int64,Int64} )
  n = length( assign )
  atmp = [ assign[i] for i in filter(i->(assign[i]!=0),collect(1:n))]
  setdiff( collect(1:n), atmp )
end
  
function perm_exact( H::Array{Float64,2} )
  n = size( H )[1]
  P = permutations(collect(1:n))
  mxscores = [ (sum( H[i,p[i]] for i = 1:n ),[p[i] for i = 1:n]) for p in P ]
  #println("mxscores: ",mxscores)
  findmaxall(mxscores)[1]
end
