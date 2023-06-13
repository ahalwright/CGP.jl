# Implementing promiscuous G-P maps via a matrix.  See the Overleaf project "Promiscuous G-P maps".
using LinearAlgebra  # enables I = identity matrix

# Returns rational matrix corresponding to Figure 2 of Garcia-Galindo paper.
function garcia_example_matrix()
  [1//2 0 1//4; 
   0 3//4 1//2; 
   1//4 0 1//4; 
   1//4 1//4 0]
end

# Temporarily says any genotype can mutate to any genotype 
function garcia_mutation_matrix()
  Rational{Int64}
  [0//1 1//1 1//1; 
   1//1 0//1 1//1; 
   1//1 1//1 0//1] 
end

# Returns rational matrix corresponding to a sampling example with K=2=alphabet size and L=2
# No promiscuity
function sampling_example_matrix()
  [1//1 0 0 0; 0 1//1 0 0; 0 0 1//1 0; 0 0 0 1//1]
end

# sampling example mutation matrix with K=2=alphabet size and L=2
function sampling_mutation_matrix( )
  Rational{Int64}
  [1 1 1 0; 
   1 1 0 1; 
   1 0 1 1; 
   0 1 1 1]
end

#= for geno corresponding to column g
function prom_geno_robustness( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, g::Int64 )
  ngenos = size(G)[2]
  #G[:,g]'*[N[g,:]'*G[i,:] for i=1:size(G)[1]]/ngenos
  sum( G[p,g]*sum( N[g,gp]*G[p,gp] for gp = 1:size(G)[2] ) for p = 1:size(G)[1] )/ngenos
end
=#

function prom_geno_robustness( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, g::Int64 )
end  

# Implementing promiscuous G-P maps via a matrix.  See the Overleaf project "Promiscuous G-P maps".

# Returns G matrix corresponding to an example with K=2=alphabet size and L=2
function G22_example_matrix()
  Rational{Int64}
  [1 0 0 0; 
   0 1 0 0; 
   0 0 1 0;
   0 0 0 1]
end

# sampling N mutation matrix with K=2=alphabet size and L=2
function N22_mutation_matrix( )
  Rational{Int64}
  [1 1 1 0; 
   1 1 0 1; 
   1 0 1 1; 
   0 1 1 1]
end

# Returns G matrix corresponding to an example with K=2=alphabet size and L=3
function G23_example_matrix()
Rational{Int64}[
1 0 0 0 0 0 0 0;
0 1 1 1 0 0 0 0;
0 0 0 0 1 1 1 0;
0 0 0 0 0 0 0 1]
end

# sampling unitation N mutation matrix with K=2=alphabet size and L=3 assuming null mutation
# subtract I to disable null mutation
function N23_mutation_matrix( )
Rational{Int64}[
1 1 1 1 0 0 0 0;
1 1 0 0 1 1 0 0;
1 0 1 0 1 0 1 0;
1 0 0 1 0 1 1 0;
0 0 1 1 1 0 0 1;
0 1 0 1 0 1 0 1;
0 1 1 0 0 0 1 1;
0 0 0 0 1 1 1 1]
end

function frequency( G::Matrix{Rational{Int64}} )
 total = sum(G)
 [ sum( G[ph,:]/total ) for ph = 1:size(G)[1] ]
end

# for geno corresponding to column g
function prom_geno_robustness( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, g::Int64 )
  ngenos = size(G)[2]
  res1 = G[:,g]'*[N[g,:]'*G[i,:] for i=1:size(G)[1]]/ngenos
  res2 = sum( G[p,g]*sum( N[g,gp]*G[p,gp] for gp = 1:size(G)[2] ) for p = 1:size(G)[1] )/ngenos
  (res1,res2)
end

function delta( x::Number, y::Number )
  x == y ? 1 : 0
end

function my_ones( n::Int64 )
  ones(Rational{Int64},n)
end

function row_dup_matrix( nrows::Int64, v::Vector{Rational{Int64}} )
  M = zeros(Rational{Int64},n,n)
  for i = 1:nrows
    M[i,:] = v
  end
  M
end


function prom_geno_evolvbability( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, g::Int64 )  # formula 6 of Garcia
  p = findfirst( x->x==1, G[:,g] )
  sum( (1-prod( N[g,gp]*(1-G[pp,gp]) for gp = 1:size(G)[2]))*(1-delta(p,pp)) for pp = 1:size(G)[1] )
end

function g_evolvability( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, g::Int64 )
  p = findfirst( x->x==1, G[:,g] )
  println("g: ",g,"  p: ",p)
  ssum = 0
  for pp = 1:size(G)[1]
    pprod = 1
    for gp = 1:size(G)[2]
      pprod *= N[g,gp]*(1-G[pp,gp])  # gp is a mut neighbor of g and G(gp) != pp
      println("gp: ",gp,"  N[g,gp]*(1-G[pp,gp]): ",N[g,gp]*(1-G[pp,gp]),"  pprod: ",pprod)
    end
    ssum += (1-pprod)*(1-delta(p,pp))  # pp != p
    println("pp: ",pp,"  (1-pprod)*(1-delta(p,pp)): ",(1-pprod)*(1-delta(p,pp)))
  end
ssum
end

function prom_pheno_evolvbability( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, p::Int64 )  # formula 7 of Garcia
  G[p,:]*sum( (1-prod(N[g,gp]*G[pp,gp] for gp = 1:size(G)[2]))*(1-delta(p,pp)) for pp = 1:size(G)[1] )
end

function my_exact_G_matrix()    # G[p,g] = P(p|g).  Thus, each column should sum to 1 (which it does)
  [1//1 1//1 0//1 0//1 0//1 0//1; 
   0//1 0//1 1//1 1//1 1//1 0//1; 
   0//1 0//1 0//1 0//1 0//1 1//1]
end
  
# Assumes that G is binary
function geno_evolvability( G::Matrix{Rational{Int64}}, N::Matrix{Rational{Int64}}, g::Int64 )  # formula 6 of Garcia
  function gdelta( gp::Int64, pp::Int64) 
    G[ pp, gp ] == 1
  end
  p = findfirst(x->x==1,G[:,g])
  sum( (1 - prod( N[g,gp]*(1-gdelta(gp,pp)) for gp = 1:size(G)[2] ))*(1-delta(p,pp)) for pp = 1:size(G)[1] )
end
 
function my_exact_N_matrix()
  [1//1 1 1 0 0 0; 
      0 1 0 1 1 0; 
      1 0 0 0 0 1; 
      0 1 0 0 1 1;
      0 1 0 0 1 0;
      0 0 1 0 1 1]
end
      
