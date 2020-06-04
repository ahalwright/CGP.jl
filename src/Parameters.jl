
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

function print_parameters( p::Parameters )
  println("MyInt: ",MyInt)
  println("numinputs: ",p.numinputs)
  println("numoutputs: ",p.numoutputs)
  println("numinteriors: ",p.numinteriors)
  println("numlevelsback: ",p.numlevelsback)
  println("nodearity: ",p.nodearity)
end

