# Evolvability using a dictionary to keep track of the phenotypes that contribute to evolution evolvability.
# Additional objectives are to simplify the logical structure form geno_complexity() and to
#   use neutral_evolution() and lambda_evolution() instead of mut_evolve()
export run_evo_dict, print_matrix, lambda_evolution

# ncircuits is the number of circuits (chromosomes) per goal used to compute evolvability
# ntries is the number of evolution attempts done in trying to evolve each chromosome
# maxsteps is the maximum number of steps for each evolution
# mutrate is the mutation rate if lambda_evolution is used.  
# mutrate<=0 means use neutral_evolution
# As of 7/3/21 does not work for cartesian==false because mutata_all() is not definined in LinCircuit.jl, and
# it looks like the definition of circuit_int() in LinCircuit.jl overrides the definition in Chromosome.jl 
#   (so I turned off the export in LinCircuit.jl).
function run_evo_dict( p::Parameters, gl::GoalList, ncircuits::Int64, numtries::Int64, maxsteps::Int64, mutrate::Float64=0.0;
      cartesian::Bool=true, csvfile::String="" )
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
  df.count_circs = Int64[]
  df.evo_count = Int64[]
  #println("size(df): ",size(df))
  results = pmap( g->evo_dict( p, g, ncircuits, numtries, maxsteps, mutrate, cartesian=cartesian ), gl )
  for res in results
    result_list = res[2:end]  # All results from evo_dict() except goal. Last element is evdict
    #println("result_list: ",result_list)
    param_list = [ res[1][1], p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, ncircuits, numtries, maxsteps, mutrate ]
    df_list = vcat( result_list[1:(end-1)], length(result_list[end]))  # length(result_list[end]) is evo_count
    #println("df_list: ",df_list)
    dfrow = vcat(param_list,df_list)
    #println("length(dfrow): ",length(dfrow))
    #println("dfrow: ",dfrow)
    push!(df,vcat(param_list,df_list))
  end 
  evdict = Dict{MyInt, Int64 }[]
  evset = Set{MyInt}[]
  evlist = Vector{Int64}[]
  pheno_matrix = zeros(Int64,length(gl),length(gl))
  gli = map( x->x[1], gl )
  println("gli: ",gli)
  for res in results
    push!(evdict,res[end])
    push!(evset,Set(keys(res[end])))
    evo_list = sort([Int64(k) for k in keys(res[end])])
    #println("evo_list: ",evo_list)
    push!(evlist,evo_list)
    goal_index = findfirst( x->x==res[1], gl )
    update_pheno_matrix( res[1], gli, res[end], pheno_matrix )
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
  print_matrix(pheno_matrix)
  phdf_hex = matrix_to_dataframe( pheno_matrix, gli, hex=true )
  phdf_dec = matrix_to_dataframe( pheno_matrix, gli, hex=false )
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile[1:end-4]*"_phmatrix_hex.csv", "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# ncircuits: ",ncircuits)
      CSV.write( f, phdf_hex, append=true, writeheader=true )
    end
    open( csvfile[1:end-4]*"_phmatrix_dec.csv", "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# ncircuits: ",ncircuits)
      CSV.write( f, phdf_dec, append=true, writeheader=true )
    end
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
  (intsdf,phdf_hex,df)
end

# Same as run_evo_dict, only for just one goal.
# if mutrate <= 0.0 use neutral_evolution, otherwise use lambda_evolution
function evo_dict( p::Parameters, g::Goal, ncircuits::Int64, numtries::Int64, maxsteps::Int64, mutrate::Float64; 
    cartesian::Bool=true )
  funcs = cartesian ? default_funcs(p.numinputs) : lin_funcs(p.numinputs)
  # evdict counts the number of occurences of the phenotypes that are mutationally adjacent to g
  evdict = Dict{MyInt, Int64 }()
  # circuit_dict is indexed by the integer representation of the circuit.
  circuit_dict = cartesian ? Dict{Int64,Int64}() : Dict{Vector{Int64},Int64}()
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
      (new_c,steps) = neutral_evolution( c, funcs, g, maxsteps )
    else
      (new_c,steps) = lambda_evolution( c, g, maxsteps, mutrate )
    end
    total_steps += steps
    if steps < maxsteps
      if haskey(circuit_dict,circuit_int(new_c))
        circuit_dict[circuit_int(new_c)] += 1
      else
        circuit_dict[circuit_int(new_c)] = 1
      end
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
  count_circs = length(keys(circuit_dict))
  return [ g, successes, t, total_steps, robust_sum, count_circs, evdict ]
end

function update_pheno_matrix( g::Goal, gli::Vector{MyInt}, evdict::Dict{MyInt,Int64}, pheno_matrix::Array{Int64,2} )
  #println("update_pheno_matrix: goal: ",g,"  gli: ",gli, "  evdict: ",evdict)
  println("update_pheno_matrix: goal: ",g,"  gli: ",gli )
  goal_index = findfirst( x->x==g[1], gli )
  for k in keys(evdict)
    ph_index = findfirst( x->x==k, gli )
    if ph_index != nothing
      #println("k: ",k,"  goal_index: ",goal_index,"  ph_index: ",ph_index)
      pheno_matrix[goal_index,ph_index] = evdict[k]
    end
  end
  #print_matrix(pheno_matrix)
end

function print_matrix( f::IO, matrix::Array{Int64,2} )
  (rows,columns) = size(matrix)
  for r = 1:rows
    for c = 1:columns
      @printf(f,"%4d ",matrix[r,c])
    end
    println(f)
  end
end

function print_matrix( matrix::Array{Int64,2} )
  print_matrix( Base.stdout, matrix )
end

# Convert matrix to a dataframe with a goal column and with column names taken from gli.
function matrix_to_dataframe( pheno_matrix::Array{Int64,2}, gli::Vector{MyInt}; hex::Bool=true)
  (rows,columns) = size(pheno_matrix)
  @assert rows==length(gli)
  df = DataFrame()
  if hex
    rnames = map(x->@sprintf("0x%x",x), gli )
  else
    rnames = map(x->@sprintf("%d",x), gli )
  end
  df.goal = rnames
  for c = 1:columns
    insertcols!(df,size(df)[2]+1,Symbol(rnames[c])=>pheno_matrix[:,c])
  end
  df
end
  
function lambda_evolution( c::Chromosome, g::Goal, maxsteps::Integer, mutrate::Float64 )
  p = c.params
  p.mutrate = mutrate
  funcs = default_funcs(p.numinputs)
  poisson_lambda = p.mutrate*p.numinteriors
  X = Poisson( poisson_lambda )
  step = 0
  ov = output_values( c)
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  #print("ov: ",ov,"  cur_dist: ",current_distance,"  "); print_circuit(c)
  while step < maxsteps && ov != g
    chrome_list = Chromosome[]
    dist_list = Float64[]
    for i = 1:p.lambda
      step += 1
      num_mutations =  rand(X)
      new_c = deepcopy(c) 
      for m = 1:num_mutations
        mutate_chromosome!( new_c, funcs )
      end
      new_ov = output_values( new_c)
      new_dist = hamming_distance( new_ov, g, c.params.numinputs )
      #print("new_ov: ",new_ov,"  new_dist: ",new_dist,"  "); print_circuit(new_c)
      push!(chrome_list,new_c)
      push!(dist_list,new_dist)
    end  
    (best_dist,ind) = findmin( dist_list )
    if best_dist <= current_distance
      c = chrome_list[ind]
      ov = output_values( c )
      current_distance = hamming_distance( ov, g, c.params.numinputs )
      #print("ov: ",ov,"  cur_dist: ",current_distance,"  "); print_circuit(c)
    end
  end # while
  if step == maxsteps
    println("lambda evolution failed with ",step," steps for goal: ",g)
    return (c, step)
  else
    println("lambda evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
end
