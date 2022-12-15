# Code based on "Enumeration and generation with a string automata representation"
#   by Marco Almeida, Nelma Moreira, Rogerio Reis, (2007).

function construct_delta()  # Example of page 95
  # Alphabet is {a,b,c} corresponding to [1,2,3]
  # Does not order the states using the transition table as described in the paper
  # Thus, states are {A,C,B,D} corresponding to [1,2,3,4]  # Ordering manually
  delta = Dict{Tuple{Int64,Int64},Int64}()  # delta maps (state,alphabet) pairs to states
  # A state is 1, B state is 2, C state is 3, D state is 4
  delta[(1,1)] = 3
  delta[(1,2)] = 2
  delta[(1,3)] = 1
  delta[(2,1)] = 4
  delta[(2,2)] = 1
  delta[(2,3)] = 2
  delta[(3,1)] = 2
  delta[(3,2)] = 4
  delta[(3,3)] = 1
  delta[(4,1)] = 3
  delta[(4,2)] = 4
  delta[(4,3)] = 2
  delta
end

# See page 95
# Assume the first state q0 = 1
function state_sequence( n::Int64, k::Int64, delta::Dict{Tuple{Int64,Int64},Int64} )
  state_seq = zeros(Int64,n)
  state_seq[1] = 1
  m = 1  # The order of states has been chosen up to m
  while m < n
    for s in state_seq[1:m]
      j = 1
      while j <= k
        d = s > 0 ? delta[(s,j)] : -1
        println("m: ",m,"  s: ",s,"  j: ",j,"  delta[(s,j)]: ", d )
        if !(delta[(s,j)] in state_seq[1:m])
          m = m + 1
          state_seq[m] = delta[(s,j)]
          println("m: ",m,"  ss: ",state_seq)
          j = k  # break out of while j < k loop
        end
        j += 1
      end
    end
    # m += 1
  end
  state_seq
end
 
function string_representation( n::Int64, k::Int64, delta::Dict{Tuple{Int64,Int64},Int64} )
  st_seq = state_sequence( n, k, delta )
  str_rep = zeros(Int64,n*k)
  i = 1
  for st = 1:n
    for j = 1:k
      println("i: ",i,"  st: ",st,"  j: ",j,"  st_seq[st]: ",st_seq[st],"  st_seq[delta[(st_seq[st],j)]]: ", st_seq[delta[(st_seq[st],j)]] )
      str_rep[ i ] = st_seq[delta[(st_seq[st],j)]]
      i += 1
    end
  end
  @assert map(s->s-1,str_rep) == [1,2,0,2,3,0,3,0,2,1,3,2]  # Convert to zero-based indexing, assert same as paper
  str_rep
end

function phi( n::Int64, k::Int64, delta::Dict{Tuple{Int64,Int64},Int64} ) 
  phi_rep = zeros( Int64, n )
  phi_rep[1] = 1
  phi_inv = zeros( Int64, n )
  phi_inv[1] = 1
  i = 1
  s = 1
  while s <= i
    for j = 1:k
      if !(delta[ (phi_inv[s],j) ] in phi_inv[1:i])
        phi_rep[ delta[ (phi_inv[s], j) ] ] = i+1
        phi_inv[i+1] = delta[ (phi_inv[s], j) ]
        i += 1
      end
      println("(s,j,i): ",(s,j,k),"  phi_rep: ",phi_rep,"  phi_inv: ", phi_inv )
    end
    s += 1
  end
  ( phi_rep, phi_inv ) 
end

# A over-simplified version of intially_connected.
# Assumes that in the definition of initially connected, the sequence [ q'[i] for i = 1:j ] is collect(1:j)
#   and the sequence [sigma[i] for i in 1:j ] is collect(1:j).
function simple_initially_connected( n::Int64, k::Int64, delta::Dict{Tuple{Int64,Int64},Int64} )
end

mutable struct DFA
  name::String
  n::Int64
  k::Int64
  states::Vector{String}
  sigma::Vector{Char}
  delta::Dict{Tuple{Int64,Int64},Int64}
  start::Int64
  final::Vector{Int64}
end

mutable struct dfaState
  name::String
end  

function dfa( ; name::String="dfa")
  dfa = DFA( name, 0, 0, Int64[], Int64[], Dict{Tuple{Int64,Int64},Int64}(), 0, [0] )
end

function setSigma!( dfa::DFA, slist::Vector{Char} )
  dfa.sigma = slist
end

function addState!( dfa::DFA, state::String )
  push!( dfa.states, state )
  dfa
end
