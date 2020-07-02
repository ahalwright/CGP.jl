export run_inf_alleles, inf_alleles, inf_alleles_to_tuple

using DataFrames
using CSV

function run_inf_alleles( nreps::Int64, 
    numinputs::IntRange, numoutputs::IntRange, numinteriors::IntRange, ngoals_rng::IntRange, levelsback::IntRange,
    popsize_rng::IntRange, max_pop_gens_rng::IntRange, tournsize::Int64, csvfile::String )
  nodearity = 2
  inf_alleles_param_list = inf_alleles_result_type[]
  df = DataFrame()
  df.numinputs = Int64[] 
  df.numoutputs = Int64[] 
  df.numints = Int64[] 
  df.levsback = Int64[] 
  #df.gl = Vector{Goal}[]
  df.ngoals = Int64[]
  df.popsize = Int64[] 
  df.tournsize = Int64[] 
  df.func_evals=Int64[]
  df.fitness = Float64[]
  df.gen_finished = Int64[] 
  iter = 1   # iteration
  for num_inputs = numinputs
    for num_outputs = numoutputs
      funcs = default_funcs(num_inputs)
      for ngoals = ngoals_rng
        gl = randgoallist(ngoals,num_inputs,num_outputs)
        for num_interiors = numinteriors
          for levsback = levelsback
            for max_pop_gens = max_pop_gens_rng 
              for popsize = popsize_rng
                for i = 1:nreps
                  p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, levsback ) 
                  #println("run inf_alleles: iter: ",iter)
                  ia = inf_alleles_result( p, gl, popsize=popsize,  max_pop_gens=max_pop_gens, tourn_size=tournsize )
                  push!(inf_alleles_param_list,ia)
                  iter += 1
                end
              end
            end
          end
        end
      end
    end
  end
  inf_alleles_result_list = pmap( x->inf_alleles(x), inf_alleles_param_list )
  #inf_alleles_result_list = map( x->inf_alleles(x), inf_alleles_param_list )
  for ia in inf_alleles_result_list
    new_row = inf_alleles_to_tuple(ia)
    Base.push!( df, new_row )
  end
  println()
  open( csvfile, "w" ) do f
    println(f,"MyInt: ",Main.CGP.MyInt)
    println(f,"funcs: ", Main.CGP.default_funcs(numinputs[end]))
    println(f,"nodearity: ",nodearity)
    CSV.write( f, df, append=true, writeheader=true )
  end                                                  
  df
end

function inf_alleles( ia::inf_alleles_result_type )
  funcs = default_funcs( ia.numinputs )
  context = construct_context(ia.numinputs)
  p = Parameters( numinputs=ia.numinputs, numoutputs=ia.numoutputs, numinteriors=ia.numints, numlevelsback=ia.levsback )
  pop = [ random_chromosome(p,funcs) for _ = 1:ia.popsize ]
  @assert ia.numinputs == p.numinputs
  for g = 1:ia.max_pop_gens
    fit_vector = zeros(Float64,ia.popsize)
    for i = 1:ia.popsize
      new_c = deepcopy(pop[i])
      mutate_chromosome!( new_c, funcs )
      output = execute_chromosome(new_c,context)
      ( new_c.fitness, matched_goals, matched_goals_list ) = goals_matched_hamming( output, ia.gl, p.numinputs)
      pop[i] = new_c
      fit_vector[i] = pop[i].fitness
    end
    #println("g:",g,"  fitv: ",fit_vector)
    maxfit = findmax( fit_vector )
    if Int(trunc(maxfit[1])) == p.numoutputs
      #println("optimum at generation ",g," = ",pop[maxfit[2]].fitness)
      ia.fitness = maxfit[1]
      ia.gen_finished = g
      ia.func_evals = g*ia.popsize
      print(".")
      #return (maxfit[1], g)
      return ia
    elseif g == ia.max_pop_gens  # return best fitness chromosome
      #println("Best fitness found: ",maxfit," after ",max_pop_gens," generations" )
      ia.fitness = maxfit[1]
      ia.gen_finished = ia.max_pop_gens
      ia.func_evals = ia.max_pop_gens*ia.popsize
      print(".")
      #return (maxfit[1], max_pop_gens)
      return ia
    end
    if ia.tourn_size == 0
      propsel!( pop, fit_vector )
    else
      tournsel!( pop, fit_vector, ia.tourn_size )
    end
  end
end

function inf_alleles_result( p::Parameters, gl::Vector{Goal} ;
    popsize::Int64=3, max_pop_gens::Int64=10, tourn_size::Int64=0 )
  inf_alleles_result_type( p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, gl, popsize, max_pop_gens,
        tourn_size, 0, 0.0, 0 )
end

function inf_alleles_to_tuple( ia::inf_alleles_result_type )
  (
    ia.numinputs ,
    ia.numoutputs ,
    ia.numints ,
    ia.levsback ,
    length(ia.gl) ,
    ia.popsize ,
    ia.tourn_size ,
    ia.func_evals ,
    ia.fitness ,
    ia.gen_finished 
  )
end
