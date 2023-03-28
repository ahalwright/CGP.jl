# Explore how the phenotype network adjacency matrix can be used to define an infinite population model.

function inf_pop_evolve( phnet_matrix::Matrix{Float64}, popvect::Vector{Float64}, fitvect::Vector{Float64}, ngens::Int64 )
  println("size(phnet_matrix): ",size(phnet_matrix),"  length(popvect): ",length(popvect),"  length(fitvect): ",length(fitvect))
  phn_matrix = deepcopy(phnet_matrix)/sum(phnet_matrix)
  pvect = popvect  # row vector
  #for i = 1:ngens
    new_pvect = (phn_matrix*Diagonal(fitvect))^ngens*pvect
    pvect = new_pvect/sum(new_pvect)
  #end
  vec(pvect)
end

# Better to using the built-in Julia type Diagonal which has the right behavior witout storing the zeros.
function diag_mat( vec::Vector )
  dim = length(vec)
  #println("dim: ",dim)
  mat = zeros(typeof(vec[1]),dim,dim)
  for i = 1:dim
    mat[i,i] = vec[i]
  end
  mat
end

function normalize( v::Union{Vector{Float64},Vector{Int64},Matrix{Float64},Matrix{Int64}} )
  v/sum(v)
end

function single_element_fitness( length::Int64, element_index::Union{Int64,MyInt},background_fitness::Float64 )
  fitvect = fill( background_fitness, length )
  if (typeof(element_index) <: MyInt)
    element_index += 1
  end
  fitvect[element_index] = 1.0
  fitvect
end
  
function single_element_fitness( p::Parameters, element_index::Union{Int64,MyInt},background_fitness::Float64 )
#function single_element_fitness( p::Parameters, element_index::MyInt,background_fitness::Float64 )
  fitvect = fill( background_fitness, 2^(2^p.numinputs))
  if (typeof(element_index) <: MyInt)
    element_index += 1
  end
  fitvect[element_index] = 1.0
  fitvect
end

function display_pheno_vect( vect::Vector, labels::Vector )
  @assert length( vect ) == length( labels )
  df = DataFrame( :phenos=>labels, :values=>vect )
end

function display_pheno_vects( vects::AbstractVector, row_labels::Vector )
  #@assert length( vects ) == length( row_labels )
  df = DataFrame( :labels=>row_labels )
  for i in 1:length(vects)
    vect = vects[i]
    insertcols!(df,size(df)[2]+1,"value$(i)"=>vect )
  end
  df
end

# mu is the mutation rate
function mutation( phnet_matrix::Matrix{Float64}, mu::Float64 )
  I + mu*normalize(phnet_matrix)
end
