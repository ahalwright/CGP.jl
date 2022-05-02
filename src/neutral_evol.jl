# An older simple version version of the function neutral_evolution().
# To run:
#   Go to the "CGP.jl/src" subdirectory of the code cloned from GitHub.
#   Start julia with "julia -L CGP.jl".  You should get the prompt:
# julia>
#   Now you want to include this file:
# julia> include("neutral_evol.jl")
#   Now you can set the parameters to 2 inputs, 1 output, 5 gates, 3 levelsback.  (Other values can be used.)
# julia> p = Parameters(2,1,5,3
#   And you need to set "funcs" to be the default gate functions:
# julia> funcs = default_funcs(p)
#   Then you can create a random circuit (genotype).
# julia> c = random_chromosome(p,funcs)
#   The print_circuit() function shows a somewhat condensed but raadable verion of the circuit:
# julia> print_circuit(c)
#   Then you set the target phenotype by:
# julia> ph = [0x0003]
#   You can run the neutral_evolution function by:
# julia> (ch,steps) = neutral_evol( c, funcs, 1000)
#   (Chromosome(Parameters(1, 4, 0.05, 0.0, 2, 1, 2, 5, 3), InputNode[InputNode(1, false, 0x0000), InputNode(2, false, 0x0000)], InteriorNode[InteriorNode(Func(&, 2, "AND"), Integer[2, 2], true, 0x000a), InteriorNode(Func(Main.CGP.Nor, 2, "NOR"), Integer[1, 3], true, 0x0001), InteriorNode(Func(|, 2, "OR"), Integer[3, 3], false, 0x0000), InteriorNode(Func(Main.CGP.Nand, 2, "NAND"), Integer[3, 3], true, 0x0005), InteriorNode(Func(Main.CGP.Nor, 2, "NOR"), Integer[6, 4], true, 0x000a)], OutputNode[OutputNode(7)], 0.0, 0.0), 6)
#   If the evolution succeeds, ch is the discovered circuit that maps to the target phenotype 0x0003.
#   You can print this circuit with:
# julia> print_circuit(ch)
#   And you can check that the output of the circuit is [0x0003]
# julia> output_values(ch)
#   1-element Vector{UInt16}:
#     0x000a

# Evolves a chromosome (cirucit) that maps to g starting with chromosome c.
# max_steps is the maximum number of evolutionary steps.
# If evolution hasn't succeeeded in max_steps, return nothing.
# Similar to mut_evolve except that this takes a single goal instead of a goal list as an argument.
function neutral_evol( c::Chromosome, g::Goal, max_steps::Integer, funcs::Vector{Func}; print_steps::Bool=false )
  #funcs = default_funcs( c.params.numinputs )
  step = 0
  ov = output_values( c) 
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  new_c = deepcopy(c)
  while step < max_steps && ov != g
    step += 1
    (new_c,active) = mutate_chromosome!( new_c, funcs )
    new_ov = output_values( new_c )
    new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    #println("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur dis: ",current_distance,"  new_dis: ",new_distance )
    if new_ov == ov 
      c = new_c
      if print_steps
        println("step: ",step," is neutral.")
      end
    elseif new_distance < current_distance
      if print_steps
        println("step: ",step,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance)
      end
      c = new_c
      ov = new_ov
      current_distance = new_distance
    else
      if print_steps
        print("step: ",step,"  new_output: ",new_ov,"  new circuit: ")
        print_circuit( new_c )
      end 
    end
  end
  if step == max_steps
    println("neutral evolution failed with ",step," steps for goal: ",g)
    return (nothing, step)
  else
    println("neutral evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end
