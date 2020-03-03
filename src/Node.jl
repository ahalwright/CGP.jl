export Integer, Node, InputNode, InteriorNode, OutputNode


#abstract Node     # Needed for julia v5
abstract type Node end    # Needed for julia v6


mutable struct InputNode <: Node
    index::Integer
    active::Bool
end

function InputNode(index::Integer)
    return InputNode(index, false)
end

mutable struct InteriorNode <: Node
    func::Func
    inputs::Vector{Integer}
    active::Bool
end

function InteriorNode(func::Func, inputs::Vector{Integer})
    return InteriorNode(func, inputs, false)
end

mutable struct OutputNode <: Node
    input::Integer
end
