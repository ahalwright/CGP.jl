# First written on 5/11/21.  Finally debugged on 5/8/21.
# Given a target phenotype, employs neutral (phenotype doesn't change) evolution to discover increeasingly robust genotypes that map to the target phenotype.
#using HypothesisTests
using Random

# Returns result_list which is a list of (circuit,rbst) pairs where rbst is the robustness() of the circuit.
# The first element of the list has maximum robustness.
function robust_evolution_parallel( c::Circuit, max_steps::Integer, nreps::Int64, funcs::Vector{Func} )
  g = output_values(c)
  rbst = robustness(c,funcs)
  print("par:  g: ",g,"  rbst: ",rbst,"   ")
  print_circuit(c,funcs)
  result_list = pmap( _->robust_evolution( c, max_steps ), collect(1:nreps) )
  #result_list = map( _->robust_evolution( c, max_steps ), collect(1:nreps) )
  #=
  for r in result_list
    print("r[2]: ",r[2],"   ")
    print_circuit(r[1],funcs)
  end
  =#
  sort!(result_list,lt=(x,y)->x[2]<y[2],rev=true)
  result_list
end

# returns a (circuit,rbst) pair where rbst is the robustness() of the circuit. 
function robust_evolution( c::Circuit, max_steps::Integer, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs=default_funcs(c.params)
  end
  g = output_values(c)
  rbst = robustness(c,funcs)
  print("step: ",0,"  g: ",g,"  rbst: ",rbst,"   ")
  print_circuit(c,funcs)
  for step = 1:max_steps 
    #(new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
    if typeof(c)==LinCircuit
      new_c = mutate_circuit!(deepcopy(c),funcs)
    else
      (new_c,active) = mutate_chromosome!(deepcopy(c),funcs)
    end
    new_rbst = robustness( new_c, funcs )
    new_ov = output_values( new_c )
    if new_ov != g || new_rbst < rbst
      continue
    end
    #println("step: ",step,"  new ov: ",output_values(new_c),"  new_rbst: ",new_rbst)
    c = new_c
    @assert output_values(c) == g
    rbst = new_rbst
    #print("step: ",step,"  ov: ",output_values(c),"  rbst: ",rbst,"   ")
    #print_circuit(c,funcs)
  end
  return (c,rbst)
end

# g is the target phenotype
function robust_evolution_pop( c::Circuit, popsize::Integer, max_steps::Integer, funcs::Vector{Func}=Func[] )
  if length(funcs) == 0
    funcs=default_funcs(c.params)
  end
  g = output_values( c )
  print("g: ",g,"  popsize: ",popsize,"  max_steps: ",max_steps," length(funcs): ",length(funcs),"    ")
  print_circuit(c,funcs)
  rbst = robustness(c,funcs)
  ov = g
  print("ov: ", ov, "  rbst: ", rbst, "   " )
  print_circuit(c,funcs)
  for step = 1:max_steps 
    pop = Tuple{Chromosome,Float64}[]
    for i = 1:popsize
      #(new_c,active) = mutate_chromosome!( deepcopy(c), funcs )
      if typeof(c)==LinCircuit
        new_c = mutate_circuit!(deepcopy(c),funcs)
      else
        (new_c,active) = mutate_chromosome!(deepcopy(c),funcs)
      end
      new_ov = output_values( new_c )
      new_rbst = robustness( new_c, funcs )
      if new_ov == ov && new_rbst >= rbst
        println("step: ",step,"  new ov: ",output_values(new_c),"  new_rbst: ",new_rbst)
        push!(pop,( deepcopy(new_c), new_rbst ))
      end
    end
    #print_pop(pop)
    if length(pop) > 1
      sort!(pop,lt=my_isless,rev=true)
    end
    if length(pop) > 0
      println("new c:  ")
      print_pop(pop)
      c = pop[1][1]
      rbst = pop[1][2]
    end
  end
  println("robust_evolution returned with phenotype: ", ov, " with robustness: ", rbst )
  @assert ov == g
  (c,rbst)
end
      
function print_pop( pop::Vector{Tuple{Chromosome,Float64}} )
  for i = 1:length(pop)
    #print( "ov: ",@sprintf("0x%04x",pop[i][2][1]),"  rbst: ",pop[i][2], "    " )
    print( "ov: ",output_values(pop[i][1]),"  rbst: ",pop[i][2], "    " )
    print_circuit(pop[i][1],funcs)
  end
end

function print_rlist( rlist )
  for r in rlist
    print("r[2]: ",r[2]," ov: ",output_values(r[1],funcs),"   ")
    print_circuit(r[1],funcs)
  end
end

# Compare x and y on the second component
function my_isless( x::Tuple, y::Tuple )
  return x[2] < y[2]
end

#=
function mutate_circuit!( c::Circuit, funcs::Vector{Func} )
  if typeof(c) == Chromosome
    mutate_chromosome!( c, funcs )
  elseif typeof(c) == LinCircuit
    mutate_circuit!( c, funcs )
  end
end
=#    

function robustness( c::Circuit, funcs::Vector{Func} )
  #print("robustness: c:  ")
  #print_circuit(c,funcs)
  c_output = output_values(c,funcs)
  outputs = mutate_all( c, funcs, output_outputs=true )
  #println("outputs[1]: ",outputs[1])
  robust_outputs = filter( x->x==c_output, outputs )
  return length(robust_outputs)/length(outputs)
end     
