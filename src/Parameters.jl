
export Parameters, default_parameters, print_parameters

struct Parameters
    mu::Integer
    lambda::Integer
    mutrate::Real
    targetfitness::Real

    numinputs::Integer
    numoutputs::Integer
    nodearity::Integer

    numinteriors::Integer
    numlevelsback::Integer
end

function Parameters(numinputs, numoutputs, numinteriors, numlevelsback)
    nodearity = 2
    mu = 1
    lambda = 4
    mutrate = 0.05
    targetfitness = 0.0
    MyInt_bits = MyIntBits( MyInt )
    #println("MyInt_bits: ",MyInt_bits)  
    if numinteriors > MyInt_bits && maxints_for_degen > MyInt_bits
      println("maxints_for_degen: ",maxints_for_degen,"  MyInt_bits: ",MyInt_bits,"  numinteriors: ",numinteriors)
      error("maxints_for_degen > MyInt_bits in function Parameters.  Run with a larger width MyInt" )
    end
    return Parameters(mu, lambda, mutrate, targetfitness, numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
end

function Parameters( ; numinputs=2, numoutputs=1, nodearity=2, numinteriors=4, numlevelsback=3 )
    mu = 1
    lambda = 4
    mutrate = 0.05
    targetfitness = 0.0
    MyInt_bits = MyIntBits( MyInt )
    #println("MyInt_bits: ",MyInt_bits)  
    if numinteriors > MyInt_bits && maxints_for_degen > MyInt_bits
      println("maxints_for_degen: ",maxints_for_degen,"  MyInt_bits: ",MyInt_bits,"  numinteriors: ",numinteriors)
      error("maxints_for_degen > MyInt_bits in function Parameters.  Run with a larger width MyInt" )
    end
    return Parameters(mu, lambda, mutrate, targetfitness, numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
end

function print_parameters(f:: IOStream, p::Parameters; comment::Bool=false )
  if !comment
    println(f,"MyInt: ",MyInt)
    println(f,"numinputs: ",p.numinputs)
    println(f,"numoutputs: ",p.numoutputs)
    println(f,"numinteriors: ",p.numinteriors)
    println(f,"numlevelsback: ",p.numlevelsback)
    println(f,"nodearity: ",p.nodearity)
  else
    println(f,"# MyInt: ",MyInt)
    println(f,"# numinputs: ",p.numinputs)
    println(f,"# numoutputs: ",p.numoutputs)
    println(f,"# numinteriors: ",p.numinteriors)
    println(f,"# numlevelsback: ",p.numlevelsback)
    println(f,"# nodearity: ",p.nodearity)
  end
end

function print_parameters(p::Parameters )
  f = Base.stdout
  print_parameters(f,p)
end

function MyIntBits( my_int::Type )
  if my_int == UInt8
    8
  elseif my_int == UInt16
    16
  elseif my_int == UInt32
    32
  elseif my_int == UInt64
    64
  elseif my_int == UInt128
    128
  else
    error("error in MyIntBits")
  end
end
