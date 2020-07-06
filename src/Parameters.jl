
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

function Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
    mu = 1
    lambda = 4
    mutrate = 0.05
    targetfitness = 0.0

    return Parameters(mu, lambda, mutrate, targetfitness, numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
end

function Parameters( ; numinputs=2, numoutputs=2, nodearity=2, numinteriors=4, numlevelsback=3 )
    mu = 1
    lambda = 4
    mutrate = 0.05
    targetfitness = 0.0
    return Parameters(mu, lambda, mutrate, targetfitness, numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
end

function print_parameters(f:: IOStream, p::Parameters )
  println(f,"MyInt: ",MyInt)
  println(f,"numinputs: ",p.numinputs)
  println(f,"numoutputs: ",p.numoutputs)
  println(f,"numinteriors: ",p.numinteriors)
  println(f,"numlevelsback: ",p.numlevelsback)
  println(f,"nodearity: ",p.nodearity)
end

function print_parameters(p::Parameters )
  f = Base.stdout
  print_parameters(f,p)
end

