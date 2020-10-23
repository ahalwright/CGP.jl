using DataFrames, Combinatorics, CSV

# maxreps is the number of iterations
# maxtries is how often to keep trying when an iteration fails.
# maxsteps is the maximum number of steps in mut_evolve.  Should be large, like =50000 for 3 inputs
function neighbor_complexity( g::Goal, p::Parameters, maxreps::Int64, maxsteps::Int64, maxtries::Int64=2*maxreps )
  funcs = default_funcs(p.numinputs)      
  chrome_complexity_sum = 0.0
  neighbor_complexity_sum = 0.0
  c_count = 0
  n_count = 0
  i = 0
  nrepeats = 0
  while i < maxtries && nrepeats < maxreps
    i += 1
    c = random_chromosome(p,funcs)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
      mut_evolve( c, [g], funcs, maxsteps )
    if step == maxsteps
      println("mut evolve failed for goal: ",goal)
      continue
    end
    #print("c: ")
    #print_build_chromosome(c)
    #println("complexity: ",complexity5(c))
    chrome_complexity_sum += complexity5(c)
    c_count += 1
    c_output = output_values(c)
    #println("goal: ",g,"  i: ",i,"  nrepeats: ",nrepeats,"  c_output: ",c_output)
    @assert sort(c_output) == sort(g)   
    output_chromes = mutate_all( c, funcs, output_chromosomes=true )
    for out_c in output_chromes
      #print("out_c: ")
      #print_build_chromosome(out_c[2])
      #println("complexity: ",complexity5(out_c[2]))
      neighbor_complexity_sum += complexity5(out_c[2])
      n_count += 1
    end
    #print("i: ",i,"  c_count: ",c_count,"  c_complex_sum: ",chrome_complexity_sum)
    #println("  n_count: ",n_count,"  n_complex_sum: ",neighbor_complexity_sum)
    nrepeats += 1
  end
  (chrome_complexity_sum/c_count, neighbor_complexity_sum/n_count)
end

function run_neighbor_complexity( ngoals::Int64, p::Parameters, maxreps::Int64, maxsteps::Int64, maxtries::Int64=2*maxreps;
    csvfile=csvfile )
  df = DataFrame()
  df.goal = Goal[]
  df.numinputs = Int64[]
  df.numoutputs = Int64[]
  df.numinteriors = Int64[]
  df.numlevelsback = Int64[]
  df.circuit_complexity = Float64[]
  df.neighbor_complexity = Float64[]
  goallist = randgoallist( ngoals, p.numinputs, p.numoutputs )
  for g in goallist
    println("g: ",g)
    (chr_complexity, nbr_complexity) = neighbor_complexity( g, p, maxreps, maxsteps, maxtries )
    push!(df, (g, chr_complexity, nbr_complexity))
  end
  df
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  open(csvfile,"w") do f
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    #println(f,"# run time in minutes: ",ttime/60)
    println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))    
    print_parameters(f,p,comment=true)
    println(f,"# maxreps: ",maxreps)
    println(f,"# maxsteps: ",maxsteps)
    println(f,"# maxtries: ",maxtries)
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end         
    
