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
function neutral_evol( c::Chromosome, funcs::Vector{Func}, g::Goal, max_steps::Integer; mut_inf_matrix::Matrix{Float64}=zeros(Float64,1,1), print_steps::Bool=false )
  default_funcs( c.params.numinputs )
  step = 0
  ov = output_values( c) 
  current_fitness = fitness_funct( p, c, [g], mut_inf_matrix=mut_inf_matrix )
  new_c = deepcopy(c)
  while step < max_steps && ov != g
    step += 1
    new_c = deepcopy(c)
    (new_c,active) = mutate_chromosome!( new_c, funcs, insert_gate_prob=0.0, delete_gate_prob=0.0 )
    new_ov = output_values( new_c )
    #new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    new_fitness = fitness_funct( p, c, [g], mut_inf_matrix=mut_inf_matrix )
    #=
    if step < 20
      print("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur fitness: ",current_fitness,"  new_fitness: ",new_fitness,"  " )
      print_circuit(new_c)
    end
    =#
    if new_ov == ov 
      c = new_c
      if print_steps
        print("step: ",step," is pheno neutral.  new_ov: ",new_ov,"  new_fitness: ",new_fitness,"  ")
        print_circuit(new_c)
      end
    elseif new_fitness == current_fitness
      c = new_c
      if print_steps
        print("step: ",step," is fitness neutral.  new_ov: ",new_ov,"  new_fitness: ",new_fitness,"  ")
        print_circuit(new_c)
      end      
    elseif new_fitness < current_fitness
      if print_steps
        println("step: ",step,"  new_output: ",new_ov," fitness improved from ",current_fitness," to ",new_fitness)
        print_circuit(new_c)
      end
      c = new_c
      ov = new_ov
      current_fitness = new_fitness
    else
      if print_steps && step <= 20
        print("step: ",step,"  new_output: ",new_ov,"  new circuit: ")
        print_circuit( new_c )
      end 
    end
  end
  if step == max_steps
    println("neutral evol failed with ",step," steps for goal: ",g)
    return (nothing, step)
  else
    println("neutral evol succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end

function accessible_path( c::Chromosome, g::Goal, max_steps::Integer, funcs::Vector{Func} )
  #  circuit_ints_list::Vector{Vector{UInt128}} )
  rand_fitvec = rand(2^numinputs)
  ov = MyInt(0)
  phfit1 = pheno_fitness(ph1,ov,numinputs,random_fit=rand_fitvec)
  phfit2 = pheno_fitness(ph2,ov,numinputs,random_fit=rand_fitvec)
  if phfit1 > phfit2
    tmp_ph = deepcopy(ph1)
    tmp_c = deepcopy( c1)
    tmp_fit = phfit1
    ph1 = ph2
    c1 = c2
    phfit1 = phfit2
    ph2 = tmp_ph
    c2 = tmp_c
    phfit2 = tmp_fit
  end
  step = 0
  ov = output_values( c1 ) 
  current_fitness = pheno_fitness( g, ov, c.params.numinputs, random_fit=rand_fitvec )
  new_c = deepcopy(c)
  while step < max_steps && ov != g
    step += 1
    (new_c,active) = mutate_chromosome!( new_c, funcs )
    new_ov = output_values( new_c )
    new = pheno_fitness( g, new_ov, c.params.numinputs, random_fit=rand_fitvec )
  end
end
  

function pheno_fitness( g::Goal, output_value::Goal, numinputs::Int64; random_fit::Vector{Float64}=Float64[])
  if length(random_fitvec) > 0
    return random_fitvec[g[1]]   # Assumes a single-component goal
  else
    #print("  hd: ",hamming_distance( output_value, g, numinputs ))
    return (1.0 - hamming_distance( output_value, g, numinputs ))
  end
end

#function find_path_between_phenotypes( ph1::Goal, c1::Circuit,  ph2:Goal, c2::Circuit, numinputs::Int64 ) end


