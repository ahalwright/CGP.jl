export Integer, Node, InputNode, InteriorNode, OutputNode, print_node_tuple, gate_tuple

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
    return InteriorNode(func, inputs, false, MyInt(0) )
end

function InteriorNode(func::Func, inputs::Vector{Integer})
    return InteriorNode(func, inputs, false, MyInt(0) )
end

mutable struct OutputNode <: Node
#struct OutputNode <: Node
    input::Integer
end

# 
function print_node_tuple( f::IO, node_vect::Vector{InputNode} )
  print(f,"(")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,node_vect[i].index,",")
  end
  if len > 1
    print(f,node_vect[len].index,")")
  else
    print(f,node_vect[len].index,",)") # add comma so 1-element tuple will be read as a tuple
  end
end

function print_node_tuple( node_vect::Vector{InputNode} )
  print_node_tuple( Base.stdout, node_vect )
end

function print_node_tuple( f::IO, node_vect::Vector{InteriorNode} )
  print(f," (")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,"(",node_vect[i].func.name,",",node_vect[i].inputs,"),")
  end
  if len > 1
    print(f,"(",node_vect[len].func.name,",",node_vect[len].inputs,"))")
  else
    # add comma so 1-element tuple will be read as a tuple
    print(f,"(",node_vect[len].func.name,",",node_vect[len].inputs,"),)")  
  end
end

function print_node_tuple( node_vect::Vector{InteriorNode} )
  print_node_tuple( Base.stdout, node_vect )
end

function print_node_tuple( f::IO, node_vect::Vector{OutputNode} )
  print(f," (")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,node_vect[i].input,",")
  end
  if len > 1
    print(f,node_vect[len].input,")")
  else
    print(f,node_vect[len].input,",)")
  end
end

function print_node_tuple( node_vect::Vector{OutputNode} )
  print_node_tuple( Base.stdout, node_vect )
end

# Outputs an vector of InteriorNodes to IO stream f in the format needed by circuit() and print_circuit() in Chromosome.jl
function gate_tuple( f::IO, node_vect::Vector{InteriorNode}, numinputs::Int64 )
  print(f," (")
  len = length(node_vect)
  for i in 1:(len-1)
    print(f,"(",numinputs+i,",",node_vect[i].func.name,",",node_vect[i].inputs[1],",",node_vect[i].inputs[2],"), ")
  end
  if len > 1
    print(f,"(",numinputs+len,",",node_vect[len].func.name,",",node_vect[len].inputs[1],",",node_vect[len].inputs[2],"))")
  else
    # add comma so 1-element tuple will be read as a tuple
    print(f,"(",numinputs+len,",",node_vect[len].func.name,",",node_vect[len].inputs[1],",",node_vect[len].inputs[2],"),)")
  end
end

#= moved to Chromosome.jl
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
  
