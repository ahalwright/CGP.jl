export Func, default_funcs

mutable struct Func
    func::Function
    arity::Integer
    #name::String
    name::AbstractString
end

Func(f::Function, a::Integer) = Func(f, a, string(f))

const AND = Func(&, 2)
const OR = Func(|, 2)
const XOR = Func(⊻,2)
const NOT = Func(~, 1)
const ZERO = Func(() -> UInt8(0), 0, "0")
const ONE = Func(() -> 0xff, 0, "1")

function default_funcs()
    return [AND, OR, NOT ]
#    return [AND, OR, NOT, ZERO, ONE]
#    return [ZERO, ONE]
#    return [AND, OR, XOR, NOT, ZERO, ONE]
end

