
function run_check_geno_exist_for_phenotypes( p::Parameters, gl::GoalList, max_evolve_tries::Int64, max_ev_steps::Int64, csvfile::String="" )
  found_list = Vector{UInt16}[]
  not_found_list = Vector{UInt16}[]
  goals_steps_list = pmap(g->check_geno_exist_for_phenotypes( p, g, max_evolve_tries, max_ev_steps ), gl )
  df = DataFrame()
  df.goal = map(x->x[1],goals_steps_list)
  df.steps = map(x->x[2],goals_steps_list)
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))                  
    print_parameters(f,p)
    println(f,"max_evolve_tries: ",max_evolve_tries)
    println(f,"max_ev_steps: ",max_ev_steps)
    println(f,not_found_list)
    CSV.write( f, df, append=true, writeheader=true )
  end
end


function check_geno_exist_for_phenotypes( p::Parameters, g::Goal, max_evolve_tries::Int64, max_ev_steps::Int64 )
  funcs = default_funcs(p.numinputs)
  i = 1
  while i < max_evolve_tries
    #println("i: ",i)
    c = random_chromosome(p,funcs)
    (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
        mut_evolve( c, [g], funcs, max_ev_steps, print_steps=false )
    if step < max_ev_steps
      break
    end
    i += 1
  end
  println("goal : ",g,"  i: ",i)
  return (g,i)
end  
    #=
    if i < max_evolve_tries
      @printf("circuit found for goal: [0x%04x]\n",g[1])
      push!(found_list,g)
    else
      @printf("circuit not found for goal: [0x%04x]\n",g[1])
      push!(not_found_list,g)
    end
    =#
