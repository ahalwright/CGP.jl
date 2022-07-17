# Defines the bitwise functions computed by logic gates.
# In current research as of 6/25/22, all gate functions used are arity 2.  Thus, NOT, IN1, IN2, ONE, and ZERO are not used as gate functions.
# VERY IMPORTANT:  Ones is set up as a global variable when default_funcs() or lin_funcs() is called.
# Common situation:  MyInt is UInt16 and numinputs = 3.  Then Ones=0x00ff is used as a mask in the computation of Nand() and Nor().
export Func, default_funcs, eval_func, lin_funcs
export Ones, AND, OR, XOR, Nand, NAND, Nor, NOR, NOT, Not, Zero, ZERO, One, ONE, In1, IN1, In2, IN2
export func_names, setup_funcs
using Main.CGP

mutable struct Func
    func::Function
    arity::Integer
    name::AbstractString
end

Func(f::Function, a::Integer) = Func(f, a, string(f))

#println("Func: p: ",p)
#println("Func: p: ",Main.CGP.p)
#println("Func:  Ones: ",Ones)
const AND = Func(&, 2, "AND")
const OR = Func(|, 2, "OR")
const XOR = Func(⊻, 2, "XOR")
Nand(x,y) = ~(x & y) & Ones   # Note that Ones is a global variable assigned when defaut_funcs() or lin_funcs() is called
const NAND = Func(Nand, 2, "NAND")
Nor(x,y) = ~(x | y) & Ones    # Note that Ones is a global variable assigned when defaut_funcs() or lin_funcs() is called  
const NOR = Func(Nor, 2, "NOR")
Not(x) = (~x) & Ones
const NOT = Func(Not, 1, "NOT")
In1(x,y) = x
const IN1 = Func(In1, 2, "IN1" )
In2(x,y) = y
const IN2 = Func(In2, 2, "IN1" )
Zero() = MyInt(0)
const ZERO = Func(Zero, 0, "0")
One() = Ones
const ONE = Func(One, 0, "1")
One() = Ones

global Ones
# Sets up Ones, ONE for the number of inputs
function setup_funcs( numinputs::Int64 )
  Ones = Main.CGP.construct_ones(numinputs)[numinputs]
  ONE = Func(One, 0, "1")
end

function default_funcs( p::Parameters )
  default_funcs( p.numinputs )
end

function default_funcs( numinputs::Int64 )
  Ones = Main.CGP.construct_ones(numinputs)[numinputs]
  ONE = Func(One, 0, "1")
  global Ones
#    return [NAND ]   # Macia's gate set
#    return [AND, XOR ]
#    return [AND, OR, XOR ]
#    return [AND, OR, NOT, ZERO, ONE]
#    return [ZERO, ONE]
#    return [AND, OR, XOR, NAND, NOR, NOT, ZERO, ONE]
#    return [AND, OR, NAND, NOR]   # Hu's gate set
    return [AND, OR, NAND, NOR, XOR]   # Raman's gate set
#    return [AND, OR, XOR, NAND, NOR, IN1, IN2 ]   # Raman's gate set plus IN1 and IN2
end

function lin_funcs( numinputs::Int64 )
  Ones = Main.CGP.construct_ones(numinputs)[numinputs]
  ONE = Func(One, 0, "1")
  global Ones
#     return [NAND ]   # Macia's gate set
#    return [AND, XOR ]
#    return [AND, OR, XOR ]
#    return [AND, OR, NOT, ZERO, ONE]
#    return [ZERO, ONE]
#    return [AND, OR, XOR, NAND, NOR, NOT, ZERO, ONE]
     return [AND, OR, NAND, NOR]   # Hu's gate set
#    return [AND, OR, XOR, NAND, NOR]   # Raman's gate set
#     return [AND, OR, XOR, NAND, NOR, IN1, IN2 ]   # Raman's gate set plus IN1 and IN2
end

# Evaluate a Func for debugging purposes
# Example:  eval_func(XOR)
function eval_func( f::Func )
    context = f.arity > 0 ? construct_contexts(f.arity)[f.arity] : []
    Main.CGP.apply( f.func, [context[i] for i = 1:f.arity ] )
end

function func_names( funcs::Vector{Func} )
  [ f.name for f in funcs ]
end
