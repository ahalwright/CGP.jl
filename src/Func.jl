export Func, default_funcs, eval_func, Ones, AND, OR, XOR, Nand, NAND, Nor, NOR, NOT, Zero, ZERO, One, ONE
using Main.CGP

mutable struct Func
    func::Function
    arity::Integer
    name::AbstractString
end

Func(f::Function, a::Integer) = Func(f, a, string(f))

#println("Func: p: ",p)
#println("Func: p: ",Main.CGP.p)
const numinputs = Main.CGP.numinputs
const Ones = Main.CGP.construct_ones(numinputs)[numinputs]
#println("Func:  Ones: ",Ones)
const AND = Func(&, 2)
const OR = Func(|, 2)
const XOR = Func(⊻,2)
Nand(x,y) = ~(x & y) & Ones
const NAND = Func(Nand,2,"NAND")
Nor(x,y) = ~(x | y) & Ones
const NOR = Func(Nor,2,"NOR")
Not(x) = (~x) & Ones
const NOT = Func(Not, 1, "NOT")
Zero() = MyInt(0)
const ZERO = Func(Zero, 0, "0")
One() = Ones
const ONE = Func(One, 0, "1")

function default_funcs()
#    return [NAND ]   # Macia's gate set
#    return [AND, OR, XOR ]
#    return [AND, OR, NOT, ZERO, ONE]
#    return [ZERO, ONE]
#    return [AND, OR, XOR, NAND, NOR, NOT, ZERO, ONE]
    return [AND, OR, XOR, NAND, NOR]   # Raman's gate set
end

# Evaluate a Func for debugging purposes
# Example:  eval_func(XOR)
function eval_func( f::Func )
    context = f.arity > 0 ? construct_contexts(f.arity)[f.arity] : []
    Main.CGP.apply( f.func, [context[i] for i = 1:f.arity ] )
end

