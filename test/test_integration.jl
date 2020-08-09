#include("CGP.jl")

# 
function integration( XI::Vector{Vector{MyInt}}, numinputs )
  X = vcat( XI... )   # Combines all of the lists in XI into a long list
  sum( entropy_lst(XI[i],numinputs) for i = 1:length(XI) ) - entropy_lst(X,numinputs)
end

function entropy_lst( X::Vector{MyInt}, numinputs )
  to_bin = to_binary( X, 2^numinputs )
  sum( entropy( to_bin[:,i]) for i = 1:(2^numinputs) )
end

# Mutual information of the subset of XI whose indices are subset_inds and the complement of this subset
function MI( XI::Vector{Vector{MyInt}}, subset_inds::Vector{Int64} )
  Xinds = collect(1:length(XI))
  entropy_lst( XI[subset_inds] ) + entropy_lst( XI[setdiff(Xinds,subset_inds)] ) - entropy_lst(XI)
end

# Mutual information of X and Y where X and Y are both lists of lists.
function MI( X::Vector{Vector{MyInt}}, Y::Vector{Vector{MyInt}} )
  entropy_lst( X ) + entropy_lst( Y ) - entropy_lst( vcat(X,Y) )
end

function complexity1( XI::Vector{Vector{MyInt}} )
  n = length(XI)
  IntXI = integration(XI)
  Xinds = collect(1:length(XI))
  ssum = 0.0
  for k = 1:(n-1)
    subsets = [XI[setdiff(Xinds,s)] for s in combinations(Xinds,k)]
    println("subsets: ",subsets)
    IXk = map(integration, subsets )
    println("k: ",k,"  IXk: ",IXk)
    ssum += (k/n)*IntXI - sum( IXk )/length(IXk)
  end
  #println("result: ",ssum)
  ssum
end

function complexity2( XI::Vector{Vector{MyInt}} )
  n = length(XI)
  Xinds = collect(1:length(XI))
  ssum = 0.0
  for k = 1:(n-1)
    subsets = [[setdiff(Xinds,s),s] for s in combinations(Xinds,k)]
    println("subsets: ",subsets)
    Xsubsets = [[XI[setdiff(Xinds,s)],XI[s]] for s in combinations(Xinds,k)]
    println("Xsubsets: ",Xsubsets)
    MIlist = [MI(Xsubsets[i][1],Xsubsets[i][2]) for i = 1:length(Xsubsets)]
    println("MIlist: ",MIlist)
  end
  sum(MIlist)/2
end

function complexity5i( )
end

function to_binary( x::MyInt, numbits::Int64 )
  result = MyInt[]
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
  result = zeros(MyInt,length(X),numbits)
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

function entropy( X::Vector{Vector{MyInt}} )
  sum( entropy(x) for x in X )
end
