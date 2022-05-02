# First written on 5/2/21.  
# Test hypothesis that choosing mutations that increase complexity will speed up evolution to a complex goal.
# Based on "Algorithmically probable mutations reproduce aspects of evolution . . . " by Santiago Hernandes-Orcozo,
#  Narsis Kiani, and Hector Zenil. http://dx.doi.org/10.1098/rsos.180399

function test_complexity_evolution( p::Parameters, gl::GoalList, popsize::Integer, max_steps::Integer, nreps::Integer, funcs::Vector{Func} )
  df = DataFrame()
  df.pheno = String[]
  df.count = Int64[]
  df.csteps_avg = Float64[]
  #df.csucceed = Bool[]
  df.dsteps_avg = Float64[]
  #df.dsucceed = Bool[]
  for g in gl
    csteps_sum = 0
    dsteps_sum = 0
    step_list = pmap(_->test_cmplx_helper( p, g, popsize, max_steps, funcs ), 1:nreps )
    #step_list = map(_->test_cmplx_helper( p, g, popsize, max_steps, funcs ), 1:nreps )
    println("step_list: ",step_list)
    csteps_sum = 0
    dsteps_sum = 0
    count = 0
    for tt in step_list
      if tt != nothing
        csteps_sum += tt[2]
        dsteps_sum += tt[3]
        count += 1
      end
    end
    println("count: ",count)
    push!(df,(@sprintf("0x%04x",g[1]),count,csteps_sum/count,dsteps_sum/count))
  end
    #=
    for r = 1:nreps
      c = random_chromosome(p,funcs)
      (cc,csteps) = complexity_evolution( deepcopy(c), g, popsize, max_steps, funcs, true )
      csucceed = cc != nothing
      (dc,dsteps) = complexity_evolution( deepcopy(c), g, popsize, max_steps, funcs, false )
      dsucceed = dc != nothing
      if csucceed && dsucceed
        count += 1
        csteps_sum += csteps
        dsteps_sum += dsteps
      end    
    end
    push!(df,(@sprintf("0x%04x",g[1]),count,csteps_sum/count,dsteps_sum/count))
  end
  =#
  df
end

function test_cmplx_helper( p::Parameters, g::Goal, popsize::Integer, max_steps::Integer, funcs::Vector{Func} ) 
  c = random_chromosome(p,funcs)
  (cc,csteps) = complexity_evolution( deepcopy(c), g, popsize, max_steps, funcs, true )
  csucceed = cc != nothing
  (dc,dsteps) = complexity_evolution( deepcopy(c), g, popsize, max_steps, funcs, false )
  dsucceed = dc != nothing
  if csucceed && dsucceed
    return (g[1],csteps,dsteps,)
  else
    return nothing
  end
end
      

function complexity_evolution( c::Chromosome, g::Goal, popsize::Integer, max_steps::Integer, funcs::Vector{Func}, use_cmplx::Bool; print_steps::Bool=false )
  my_isless =  use_cmplx ? cmplx_isless : dist_isless
  step = 0
  ov = output_values( c )
  #println("    ov: ",output_values(c))
  current_distance = hamming_distance( ov, g, c.params.numinputs ) 
  #println("cur ov: ",output_values(c),"  cur_dist: ",current_distance)
  while step < max_steps 
    pop = Tuple{Chromosome,Goal,Float64,Float64}[]
    for i = 1:popsize
      (new_c,active) = mutate_chromosome!( c, funcs )
      new_ov = output_values( new_c )
      new_distance = hamming_distance( new_ov, g, c.params.numinputs )  
      #println("new ov: ",output_values(new_c),"  new_dist: ",new_distance)
      new_cmplx = complexity5(new_c)
      push!(pop,(deepcopy(new_c),new_ov,new_distance,new_cmplx))
    end
    #print_pop(pop)
    sort!(pop,lt=my_isless)
    #println("====")
    #print_pop(pop)
    if pop[1][3] == 0.0
      println("complexity evolution succeeded with ",step," steps for goal: ",g)
      #print_circuit(pop[1][1])
      #println("g: ",g,"  output_values(pop[1][1]): ",output_values(pop[1][1]))
      @assert output_values(pop[1][1]) == g
      return (pop[1][1],step)
    elseif pop[1][3] <= current_distance
      #println("new c")
      c = pop[1][1]
    end
    step += 1
  end
  println("complexity evolution failed with ",step," steps for goal: ",g)
  return (nothing, step) 
end
      
function print_pop( pop::Vector{Tuple{Chromosome,Goal,Float64,Float64}} )
  for i = 1:length(pop)
    print( "ov: ",@sprintf("0x%04x",pop[i][2][1]),"  dis: ",pop[i][3],"  cmplx: ",pop[i][4], "    " )
    print_circuit(pop[i][1])
  end
end

# Compare x and y on the first by decreasing distance 
function dist_isless( x::Tuple, y::Tuple )
  return x[3] < y[3]
end
    
# Compare x and y on the first by decreasing distance and then by increasing complexity
function cmplx_isless( x::Tuple, y::Tuple )
  if x[3] < y[3]
    return true
  elseif x[3] == y[3] && x[4] > y[4]   # We want minimum distance and maximum complexity
    return true
  else
    return false
  end
end
