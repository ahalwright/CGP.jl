# First written on 5/2/21.  Finally debugged on 5/8/21.
# Test hypothesis that choosing mutations that increase complexity will speed up evolution to a complex goal.
# Based on "Algorithmically probable mutations reproduce aspects of evolution . . . " by Santiago Hernandes-Orcozo,
#  Narsis Kiani, and Hector Zenil. http://dx.doi.org/10.1098/rsos.180399
using HypothesisTests
using Random

function test_complexity_evolution( p::Parameters, gl_counts::Vector{Tuple{Vector{MyInt},Int64}},  popsize::Integer, max_steps::Integer, nreps::Integer, funcs::Vector{Func}; csvfile::String="" )
  df = DataFrame()
  df.pheno = String[]
  df.count = Int64[]
  df.nreps = Int64[]
  df.csteps_mean = Float64[]
  df.csteps_std = Float64[]
  df.dsteps_mean = Float64[]
  df.dsteps_std = Float64[]
  df.pvalue = Float64[]
  i = 1
  for gc in gl_counts
    println("gc: ",gc)
    count = gc[2]
    i += 1
    fslist = pmap(i->test_cmplx_helper( p, gc[1], popsize, max_steps ), 1:nreps )
    #fslist = map(i->test_cmplx_helper( p, gc[1], popsize, max_steps ), 1:nreps )
    #println("step_list: ",step_list)
    #fslist = filter( x->x!=nothing, step_list )
    #fslist = step_list 
    #println("fslist: ",fslist)
    if length(fslist) > 0
      csteps_mean = mean( map(x->x[2], fslist ))
      dsteps_mean = mean( map(x->x[3], fslist ))
      csteps_std = std( map(x->x[2], fslist ))
      dsteps_std = std( map(x->x[3], fslist )) 
      #println("csteps_mean: ",csteps_mean)
      #println("dsteps_mean: ",dsteps_mean)
      #println("csteps_std: ",csteps_std)
      #println("dsteps_std: ",dsteps_std)
      pval = 100.0
      try
        pval = length(fslist)>=4 ? pvalue( UnequalVarianceTTest( map(x->x[2],fslist), map(x->x[3],fslist))) : -1.0 
      catch
        pval = 100.0
      end
      println("pval: ",pval)
      push!(df,(@sprintf("0x%04x",gc[1][1]),gc[2],length(fslist),csteps_mean,csteps_std,dsteps_mean,dsteps_std,pval))
    else
      push!(df,(@sprintf("0x%04x",gc[1][1]),gc[2],length(fslist),nreps,0,0,0.0))
    end
  end
  if length(csvfile) > 0
    write_df_to_csv( df, p, funcs, csvfile, popsize=popsize, max_steps=max_steps, nreps=nreps )
  end
  df
end

function test_cmplx_helper( p::Parameters, g::Goal, popsize::Integer, max_steps::Integer ) 
  funcs = default_funcs(p)   # Necessary to be done on each process when called by pmap()
  c = random_chromosome(p,funcs)
  (cc,csteps) = complexity_evolution( deepcopy(c), g, popsize, max_steps, funcs, true )
  csucceed = cc != nothing
  (dc,dsteps) = complexity_evolution( deepcopy(c), g, popsize, max_steps, funcs, false )
  dsucceed = dc != nothing
  if csucceed && dsucceed
    return (g[1],csteps,dsteps,)
  elseif csucceed && !dsucceed
    return (g[1],csteps,max_steps)
  elseif !csucceed && dsucceed
    return (g[1],max_steps,dsteps,)
  else
    return (g[1],max_steps,max_steps)
  end
end

function complexity_evolution( c::Chromosome, g::Goal, popsize::Integer, max_steps::Integer, funcs::Vector{Func}, use_cmplx::Bool; print_steps::Bool=false )
  #println("g: ",g,"  popsize: ",popsize,"  max_steps: ",max_steps," length(funcs): ",length(funcs),"  use_cmplx: ",use_cmplx)
  #print_circuit(c)
  my_isless =  use_cmplx ? cmplx_isless : dist_isless
  step = 0
  ov = output_values( c )
  #println("    ov: ",output_values(c))
  current_distance = hamming_distance( ov, g, c.params.numinputs ) 
  #println("cur ov: ",ov,"  cur_dist: ",current_distance)
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
