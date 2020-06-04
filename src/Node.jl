export Integer, Node, InputNode, InteriorNode, OutputNode

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
