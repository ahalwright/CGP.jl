export Integer, Node, InputNode, InteriorNode, OutputNode, print_node_tuple

abstract type Node end  
mutable struct Func
    func::Function
    arity::Integer
    name::AbstractString
end
Func(f::Function, a::Integer) = Func(f, a, string(f)) 


mutable struct InputNode <: Node
    index::Integer
    active::Bool
    cache::MyInt
end

function InputNode(index::Integer)
    return InputNode(index, false, 0 )
end

mutable struct InteriorNode <: Node
    func::Func
    inputs::Vector{Integer}
    active::Bool
    cache::MyInt
end

function InteriorNode(func::Func, inputs::Vector{Int64})
    return InteriorNode(func, inputs, false, 0 )
end

function InteriorNode(func::Func, inputs::Vector{Integer})
    return InteriorNode(func, inputs, false, 0 )
end

mutable struct OutputNode <: Node
    input::Integer
end

function print_node_tuple( f::IO, node_vect::Vector{InputNode} )
  print(f,"  (")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,node_vect[i].index,",")
  end
  print(f,node_vect[len].index,")")
end

function print_node_tuple( node_vect::Vector{InputNode} )
  print_node_tuple( Base.stdout, node_vect )
end

function print_node_tuple( f::IO, node_vect::Vector{InteriorNode} )
  print(f,"  (")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,"(",node_vect[i].func.name,",",node_vect[i].inputs,"),")
  end
  print(f,"(",node_vect[len].func.name,",",node_vect[len].inputs,"))")
end

function print_node_tuple( node_vect::Vector{InteriorNode} )
  print_node_tuple( Base.stdout, node_vect )
end

function print_node_tuple( f::IO, node_vect::Vector{OutputNode} )
  print(f,"  (")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,node_vect[i].input,",")
  end
  print(f,node_vect[len].input,")")
end

function print_node_tuple( node_vect::Vector{OutputNode} )
  print_node_tuple( Base.stdout, node_vect )
end

#=
function print_build_chromosome( f::IO, c::Chromosome )
  println(f, "build_chromosome(")
  print_node_tuple(f, c.inputs )
  print(f,",")
  print_node_tuple(f, c.interiors )
  print(f,",")
  print_node_tuple(f, c.outputs )
  println(f, ")")
end
=#
  
