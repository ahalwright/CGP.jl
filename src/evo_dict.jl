# Evolvability using a dictionary to keep track of the phenotypes that contribute to evolution evolvability.
# Additional objectives are to simplify the logical structure form geno_complexity() and to
#   use neutral_evolution() and lambda_evolution() instead of mut_evolve()

# ncircuits is the number of circuits (chromosomes) per goal used to compute evolvability
# ntries is the number of evolution attempts done in trying to evolve each chromosome
# maxsteps is the maximum number of steps for each evolution
# mutrate is the mutation rate if lambba_evolution is used.  
# mutrate<=0 means use neutral_evolution
function run_evo_dict( p::Parameters, gl::GoalList, ncircuits::Int64, numtries::Int64, maxsteps::Int64, mutrate::Float64=0.0;
      cartesian::Bool=true, csvfile::String=""  )
  df = DataFrame()
  #df.goal = Goal[]
  df.goal = MyInt[]
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numgates = Int64[]
  df.numlevsback = Int64[]
  df.ncircuits = Int64[]
  df.numtries = Int64[]
  df.maxsteps = Int64[]
  df.mutrate = Float64[]
  df.successes = Int64[]
  df.tries = Int64[]
  df.totalsteps = Int64[]
  df.robust_sum = Int64[]  
  df.evo_count = Int64[]
  results = pmap( g->evo_dict( p, g, ncircuits, numtries, maxsteps, mutrate, cartesian=cartesian ), gl )
  for res in results
    result_list = res[2:end]
    #println("result_list: ",result_list)
    param_list = [ res[1][1], p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, ncircuits, numtries, maxsteps, mutrate ]
    #println("length(param_list): ",length(param_list))
    #println("param_list: ",param_list)
    df_list = vcat( result_list[1:(end-1)], length(result_list[end]))
    #println("df_list: ",df_list)
    push!(df,vcat(param_list,df_list))
  end 
  evdict = Dict{MyInt, Int64 }[]
  evset = Set{MyInt}[]
  evlist = Vector{Int64}[]
  for res in results
    push!(evdict,res[end])
    push!(evset,Set(keys(res[end])))
    evo_list = sort([Int64(k) for k in keys(res[end])])
    #println("evo_list: ",evo_list)
    push!(evlist,evo_list)
  end
  #println("evlist: ",evlist)
  insertcols!(df,size(df)[2]+1,:evo_phenos=>evlist)
  intsdf = DataFrame()
  for g in gl
    insertcols!(intsdf,size(intsdf)[2]+1,Symbol(@sprintf("0x00%x",g[1]))=>Int64[])
  end
  for evs1 in evset
    intersects = [ length(intersect(evs1,evs2)) for evs2 in evset ]
    println("intersects: ",intersects)
    push!(intsdf,intersects)
  end  
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile[1:end-4]*"_ints.csv", "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      println(f,"# numinputs: ",p.numinputs)
      println(f,"# numoutputs: ",p.numoutputs)
      println(f,"# numinteriors: ",p.numinteriors)
      println(f,"# numlevelsback: ",p.numlevelsback)
      println(f,"# ncircuits: ",ncircuits)
      println(f,"# numtries: ",numtries)
      println(f,"# maxsteps: ",maxsteps)
      CSV.write(f, intsdf, append=true, writeheader=true )
    end
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      println(f,"# ncircuits: ",ncircuits)
      println(f,"# numtries: ",numtries)
      println(f,"# maxsteps: ",maxsteps)
      CSV.write(f, df, append=true, writeheader=true )
    end
  end
  (intsdf,df)
end

# Same as run_evo_dict, only for just one goal.
# if mutrate <= 0.0 use neutral_evolution, otherwise use lambda_evolution
function evo_dict( p::Parameters, g::Goal, ncircuits::Int64, numtries::Int64, maxsteps::Int64, mutrate::Float64; 
    cartesian::Bool=true )
  funcs = cartesian ? default_funcs(p.numinputs) : lin_funcs(p.numinputs)
  evdict = Dict{MyInt, Int64 }()
  robust_sum = 0
  total_steps = 0
  successes = 0
  t = 0
  while t < numtries && successes < ncircuits
    t += 1
    if cartesian 
      c = random_chromosome( p, funcs )
    else
      c = rand_lcircuit(p,funcs)  
    end
    if mutrate <= 0.0
      (new_c,steps) = neutral_evolution( c, g, maxsteps )
    else
      (new_c,steps) = lambda_evolution( c, g, maxsteps, mutrate )
    end
    total_steps += steps
    if steps < maxsteps
      successes += 1
      ov = output_values( new_c )   
      @assert g == ov
      #println("g: ",g,"  ov: ",ov)
      neighboring_phenos = mutate_all( new_c, funcs, output_outputs=:true )
      #println("neighboring_phenos: ",neighboring_phenos[1:5])
      robust_outputs = filter( x->x==ov, neighboring_phenos )
      robust_sum += length(robust_outputs)
      for ph in neighboring_phenos
        if haskey(evdict,ph[1])
          evdict[ph[1]] += 1
        else
          evdict[ph[1]] = 1
        end
      end
    end  # if
  end  # while
  return [ g, successes, t, total_steps, robust_sum, evdict ]
end
