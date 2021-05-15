using DataFrames
# Test evolve_function.jl in Evolve.jl for the specific purpose of reducing the number of gates
#    in a 2-component goal where circuits has already been evolved to map to each single-component goal.

# Goal is the 2-component goal
function test_evolve_function( goal::Goal, p::Parameters, maxsteps::Int64 )
  @assert p.numouputs==2
  p1 = Parameters( div(p.numinputs,2), 1, p.numlevelsback, p.numinteriors, p.numlevelsback )
  maxtries = 10
  tries = 1
  done = :false
  while !done && tries < maxtries
    ec = random_chromosome(p1)
    (c1,steps,maxsteps) = neutral_evolution( ec, [goal[1]], maxsteps )
    if steps < maxsteps
      done = :true
    end
  end
  if !done
    error("unable to evolve circuit to map to goal ",[goal[1]])
  end
  tries = 1
  done = :false
  while !done && tries < maxtries
    ec = random_chromosome(p1)
    (c2,steps,maxsteps) = neutral_evolution( ec, [goal[2]], maxsteps )
    if steps < maxsteps
      done = :true
    end
  end
  if !done
    error("unable to evolve circuit to map to goal ",[goal[1]])
  end
end  

# Combines two 1-output circuits c1 and c2 with the same paramterers and number of gates into a 2-output circuit
# The output values of the result is the 2-output goal whose components are the output values of c1 and of c2.
function combine_circuits( c1::Chromosome, c2::Chromosome )
  # Two functions used only within this function
  inp1(i,numinputs) = i <= numinputs ? i : 2*(i-numinputs)+numinputs-1
  inp2(i,numinputs) = i <= numinputs ? i : 2*(i-numinputs)+numinputs
  @assert c1.params.numinputs == c2.params.numinputs
  @assert c1.params.numoutputs == 1
  @assert c2.params.numoutputs == 1
  if length(c1.interiors) < length(c2.interiors)
    println("  c1ints: ",length(c1.interiors),"  c2ints: ",length(c2.interiors))
    ctemp = c1
    c1 =c2
    c2 = ctemp
  end
  #println("c1ints: ",length(c1.interiors),"  c2ints: ",length(c2.interiors))
  p = Parameters( c1.params.numinputs, 2, c1.params.numinteriors+c2.params.numinteriors,
    c1.params.numlevelsback+c2.params.numlevelsback )
  inputs = c1.inputs
  interiors = InteriorNode[]
  new_numinteriors = 2*length(c2.interiors) 
  #println("new_numints: ",new_numinteriors)
  #print_parameters(p)
  for i = 1:length(c2.interiors)
    int_node1 = deepcopy(c1.interiors[i])
    input1 = inp1(c1.interiors[i].inputs[1],p.numinputs)
    input2 = inp1(c1.interiors[i].inputs[2],p.numinputs)
    int_node1.inputs = [input1,input2]
    #println("i: ",i,"  int_node1: ",int_node1)
    push!(interiors,int_node1)
    int_node2 = deepcopy(c2.interiors[i])
    input1 = inp2(c2.interiors[i].inputs[1],p.numinputs)
    input2 = inp2(c2.interiors[i].inputs[2],p.numinputs)
    int_node2.inputs = [input1,input2]
    #println("i: ",i,"  int_node2: ",int_node2)
    push!(interiors,int_node2)
  end
  #println("interiors: ",interiors)
  outputs = [OutputNode(p.numinputs+new_numinteriors-1),OutputNode(p.numinputs+new_numinteriors)]
  Chromosome( p, c1.inputs, interiors, outputs, 0.0, 0.0 )
end

function test_combine_circuits( p::Parameters )
  c1 = random_chromosome(p); print(output_values(c1),"  "); print_circuit(c1) 
  c2 = random_chromosome(p); print(output_values(c2),"  "); print_circuit(c2) 
  c = combine_circuits(c1,c2); print(output_values(c),"  "); print_circuit(c) 
end

# Evaluate the performance of function evolve_reduce_numactive() which is in evolve_funct.jl.
# For each pair of chromosomes
#   Creates two chromosomes using parameter settings p.
#   Combines them into one chromosome with one more input
#   For each trial run evolve_reduce_numactive() and record average and minimum numactive.
function test_reduce_numactive( p::Parameters, numpairs::Int64, numtrials::Int64, maxsteps::Int64; 
      num_mutations::Int64=1, csvfile::String="" )
  df = DataFrame()
  df.combined_goal = Goal[]
  df.numinputs = Int64[]
  df.numgates = Int64[]
  df.numlevelsback = Int64[]
  df.numpairs = Int64[]
  df.numtrials = Int64[]
  df.maxsteps = Int64[]
  df.nactive0 = Float64[]
  df.nactive1_mean = Float64[]
  df.nactive1_std = Float64[]
  if num_mutations >= 2
    df.nactive2_mean = Float64[]
    df.nactive2_std = Float64[]
  end
  if num_mutations == 3
    df.nactive3_mean = Float64[]
    df.nactive3_std = Float64[]
  end
  for npair = 1:numpairs
    c1 = random_chromosome(p); print(output_values(c1),"  "); print_circuit(c1) 
    c2 = random_chromosome(p); print(output_values(c2),"  "); print_circuit(c2) 
    c = combine_circuits(c1,c2)
    rowlist = [ output_values(c), p.numinputs, p.numinteriors, p.numlevelsback, numpairs, numtrials, maxsteps ]
    #@assert [output_values(c1)[1],output_values(c1)[1]] == output_values(c)
    println("output_values(c): ",output_values(c))
    println("goal: ",[output_values(c1)[1],output_values(c2)[1]])
    #nactive_lists = map(i->reduce_nactive_helper(c,maxsteps,num_mutations=num_mutations), collect(1:numtrials ) )
    nactive_lists = pmap(i->reduce_nactive_helper(c,maxsteps,num_mutations=num_mutations), collect(1:numtrials ) )
    nactive0 = vcat( [ na[1] for na in nactive_lists]... ) 
    nactive1 = vcat( [ na[2] for na in nactive_lists]... )
    if num_mutations == 1
      row_results = [mean(nactive0),mean(nactive1),std(nactive1)]
    end
    if num_mutations >= 2
      nactive2 = vcat( [ na[3] for na in nactive_lists]... )
      row_results = [mean(nactive0),mean(nactive1),std(nactive1),mean(nactive2),std(nactive2)]
    end
    if num_mutations >= 3
      nactive3 = vcat( [ na[4] for na in nactive_lists]... )
      row_results = [mean(nactive0),mean(nactive1),std(nactive1),mean(nactive2),std(nactive2),mean(nactive3),std(nactive3)]
    end
    row = deepcopy(rowlist)
    println("length(row): ",length(row),"  length(row_results): ",length(row_results),"  size(df): ",size(df))
    #println("npair: ",npair,"  mean(nactive1): ",mean(nactive1),"  mean(nactive2): ",mean(nactive2),"  mean(nactive3): ",mean(nactive3))
    push!(df,vcat(row,row_results))
  end
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# numpairs: ", numpairs)
      println(f,"# numtrials: ", numtrials)
      println(f,"# maxsteps: ", maxsteps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end 

# nactive0 is number_active_gates() after evolution to find circuit that maps to goal
#    but before running evolve_reduce_numactive()
function reduce_nactive_helper( c::Chromosome, maxsteps::Int64; num_mutations::Int64=1 ) 
  #println("num_mutations: ",num_mutations)
  goal = output_values( c )
  (ec,steps) = neutral_evolution( random_chromosome(c.params), goal, maxsteps )
  nactive0 = number_active_gates(ec)
  println("nactive0: ",nactive0)
  nactive1 = Int64[]
  nc = deepcopy(c)
  (ec,cur_nactive) = evolve_reduce_numactive( nc, maxsteps, num_mutations=1 ) 
  @assert cur_nactive == number_active_gates(ec)
  push!(nactive1,number_active_gates(ec))
  if num_mutations == 1
    return(nactive0, nactive1)
  end
  nactive2 = Int64[]
  nc = deepcopy(c)
  (ec,cur_nactive) = evolve_reduce_numactive( nc, maxsteps, num_mutations=2 ) 
  push!(nactive2,number_active_gates(ec))
  if num_mutations == 2
    return(nactive0, nactive1, nactive2)
  end
  nactive3 = Int64[]
  nc = deepcopy(c)
  (ec,cur_nactive) = evolve_reduce_numactive( nc, maxsteps, num_mutations=3 ) 
  push!(nactive3,number_active_gates(ec))
  #println("trial: ",trial,"  mean(nactive1): ",mean(nactive1),"  mean(nactive2): ",mean(nactive2),"  mean(nactive2): ",mean(nactive2))
  (nactive0, nactive1, nactive2, nactive3)
end
  
# folder_path is relative to data folder 
function compare_dataframe( folder_path::String, suffix::String )
  rdf = read_dataframe("$(folder_path)/reduce_numactive_$(folder_path)_$(suffix).csv")
  kdf = read_dataframe("$(folder_path)/k_complexity$(folder_path)$(suffix).csv")
  cdf = DataFrame()
  cdf.goal = kdf.goal
  cdf.kng = kdf.num_gates
  cdf.rng = rdf.nactive1_mean
  cdf.k_minus_r_diff = cdf.kng .- cdf.rng
  cdf
end
