# Code based on "Enumeration and generation with a string automata representation"
#   by Marco Almeida, Nelma Moreira, Rogerio Reis, (2007).

mutable struct DFA
  name::String
  n::Int64
  k::Int64
  states::Vector{String}
  sigma::Vector{Char}
  delta::Dict{Tuple{Int64,Char},Int64}
  start::Int64
  final::Vector{Int64}
end

#=
function state_sequence( ddfa::DFA )
  n = ddfa.n
  k = ddfa.k
  state_seq = zeros(Int64,n)
  state_seq[1] = ddfa.start
  m = 1
  while m < n
    for s in state_seq[1:m]
      j = 1
      while j <= k
        d = s > 0 ? ddfa.delta[(s,ddfa.sigma[j])] : -1
        println("m: ",m,"  s: ",s,"  j: ",ddfa.sigma[j],"  ddfa.delta[(s,j)]: ", d )
        if !(ddfa.delta[(s,ddfa.sigma[j])] in state_seq[1:m])
          m = m + 1
          state_seq[m] = ddfa.delta[(s,ddfa.sigma[j])]
          println("m: ",m,"  ss: ",state_seq)
          j = k  # break out of while j < k loop
        end
        j += 1 
      end
      if m < n
        m += 1
      else
        break
      end
    end
  end
  state_seq
end
=#

# See page 95 of Almeida et al. 2007
# Based on Wikipedia depth first search
function state_sequence( ddfa::DFA )
  n = ddfa.n
  k = ddfa.k
  st_seq = Int64[]
  explored = fill(false,n)
  queue = Queue{Int}()
  explored[ddfa.start] = true
  enqueue!(queue,ddfa.start)
  push!(st_seq,ddfa.start)
  while length(queue) > 0
    v = dequeue!(queue)
    for j = 1:k
      w = ddfa.delta[(v,ddfa.sigma[j])]
      if !explored[w]
        explored[w] = true
        push!(st_seq, w )
        enqueue!(queue,w)
      end
    end
  end
  st_seq 
end
 
# Checks that the string representation of ddfa is the same as the example of page 95 of Almeida et al. 2007
function string_representation( ddfa::DFA )
  st_seq = state_sequence( ddfa )
  str_rep = zeros(Int64,ddfa.n*ddfa.k)
  i = 1
  for st = 1:ddfa.n
    for j in ddfa.sigma
      println("i: ",i,"  st: ",st,"  j: ",j,"  st_seq[st]: ",st_seq[st],"  st_seq[ddfa.delta[(st_seq[st],j)]]: ", st_seq[ddfa.delta[(st_seq[st],j)]] )
      str_rep[ i ] = st_seq[ddfa.delta[(st_seq[st],j)]]
      i += 1
    end
  end
  @assert map(s->s-1,str_rep) == [1,2,0,2,3,0,3,0,2,1,3,2]  # Convert to zero-based indexing, assert same as paper
  str_rep
end

function phi( ddfa::DFA) 
  phi_rep = zeros( Int64, ddfa.n )
  phi_rep[1] = 1
  phi_inv = zeros( Int64, ddfa.n )
  phi_inv[1] = 1
  i = 1
  s = 1
  while s <= i
    for j in ddfa.sigma
      if !(ddfa.delta[ (phi_inv[s],j) ] in phi_inv[1:i])
        phi_rep[ ddfa.delta[ (phi_inv[s], j) ] ] = i+1
        phi_inv[i+1] = ddfa.delta[ (phi_inv[s], j) ]
        i += 1
      end
      println("(s,j,i): ",(s,j,i),"  phi_rep: ",phi_rep,"  phi_inv: ", phi_inv )
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

function dfa( ; name::String="dfa")
  dfa = DFA( name, 0, 0, Int64[], Int64[], Dict{Tuple{Int64,Char},Int64}(), 0, [0] )
end

function setSigma!( dfa::DFA, slist::Vector{Char} )
  dfa.sigma = slist
  dfa.k = length(slist)
end

function addState!( dfa::DFA, state::String )
  push!( dfa.states, state )
  dfa.n += 1
end

function setInitial!( dfa::DFA, init::Int64 )
  dfa.start = init
end

function setFinal!( dfa::DFA, final::Vector{Int64} )
  dfa.final = final
end

function evalWordP( dfa::DFA, input::String )
  st = dfa.start
  for c in input
    st = dfa.delta[( st, c )] 
  end
  #st in dfa.final
  st
end

function addTransition!( dfa::DFA, fromSt::Int64, input::Char, toSt::Int64 )
  dfa.delta[(fromSt,input)] = toSt
end

# For example of Almeida 2007 p 95
function construct_ddfa()
  ddfa = dfa( name="ddfa" )
  setSigma!( ddfa, ['0', '1', '2'] )
  ddfa.k = length(ddfa.sigma)
  addState!( ddfa, "S1" )
  addState!( ddfa, "S2" )
  addState!( ddfa, "S3" )
  addState!( ddfa, "S4" )
  ddfa.n = length(ddfa.states)
  setInitial!(ddfa,1)
  setFinal!(ddfa,[1])
  # '0' corresponds to 'a', '1' corresponds to 'b', '2' corresponds to 'c',
  addTransition!( ddfa, 1, '0', 3)
  addTransition!( ddfa, 1, '1', 2)
  addTransition!( ddfa, 1, '2', 1)
  addTransition!( ddfa, 2, '0', 4)
  addTransition!( ddfa, 2, '1', 1)
  addTransition!( ddfa, 2, '2', 2)
  addTransition!( ddfa, 3, '0', 2)
  addTransition!( ddfa, 3, '1', 4)
  addTransition!( ddfa, 3, '2', 1)
  addTransition!( ddfa, 4, '0', 3)
  addTransition!( ddfa, 4, '1', 4)
  addTransition!( ddfa, 4, '2', 2)
  ddfa
end 

function test_DFA()
  ddfa = construct_ddfa()
  println("evalWordP(ddfa,\"011\"): ",evalWordP(ddfa,"011"))
  println("evalWordP(ddfa,\"1011\"): ",evalWordP(ddfa,"1011"))
  ddfa
end
