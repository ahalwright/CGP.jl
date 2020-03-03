export Parameters, default_parameters

struct Parameters
    mu::Integer
    lambda::Integer
    mutrate::Real
    targetfitness::Real
    mask::Integer

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
    mask = 0x1
    for i = 0:nodearity-1
      mask = mask | (mask << 2^i)
    end

    return Parameters(mu, lambda, mutrate, targetfitness, mask, numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
end

