# Given a genotype ch, does up to max_iterations where each iteration replaces ch with a higher fitness mutation
# Iterations continue until no higher fitness mutation is found or max_iterations are done.
# Returns the highest fitenss genotype found.
function hill_climb( fitness::Vector{Float64}, ch::Chromosome, funcs::Vector{Func}, max_iterations::Int64)::Tuple{Chromosome,Float64,Goal}
  fitfunct( g::Goal ) = fitness[ g[1]+1 ]::Float64   # fitness of a phenotype
  fitfunct( ch::Circuit ) = fitfunct( output_values(ch) )::Float64  # fitness of a genotype
  local current_ch
  best_ch = ch
  best_fitness = fitfunct(ch)
  # println("best_fitness: ", best_fitness)
  for i = 1:max_iterations   # Note return exit from loop
    #println("i: ",i,"  best_fitness: ",best_fitness)
    # Find the most fit out of best_ch and the mutations of best_ch 
    mutations = mutate_all(best_ch,funcs,output_circuits=true,output_outputs=false)
    push!(mutations,best_ch)
    (current_fitness,index) = findmax(map(ch->fitfunct(ch),mutations))
    current_ch = mutations[index]
    #println("i: ",i,"  output_values(current_ch): ",output_values(ch))
    #println("i: ",i,"  current_fitness: ",current_fitness)
    # Test if current_ch is better than ch.
    if best_fitness < current_fitness  # if so, replace best_ch by current_ch and continue loop
      best_ch = current_ch 
      best_fitness = fitfunct(best_ch)
      #println("i: ",i,"  new best_fitness: ", best_fitness)
    else    # if not, the hill climb is done
      #println("break:  current_fitness: ",current_fitness,"  best_fitness: ",best_fitness)
      break
    end
  end
  return (current_ch, fitfunct(current_ch),output_values(current_ch))
end

function run_hill_climb( fitness::Vector{Float64}, ch::Chromosome, funcs::Vector{Func}, max_iterations::Int64, nreps::Int64; csvfile::String="" )
  df = DataFrame( iteration=Int64[], goal=Vector{MyInt}[], fitvalue=Float64[] )
  fit_goal_list = Tuple{Float64,Goal}[]
  for i = 1:nreps
    ( fitvalue, goal ) = hill_climb( fitness, random_chromosome(ch.params,funcs), funcs, 10 )[2:3]
    row = (i,goal,fitvalue)
    push!(df,row)
    #push!( fit_goal_list, ( fitvalue, goal ) )
  end
  if length(csvfile) > 0
    println("csvfile: ",csvfile)
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs )
      print_parameters(f,p,comment=true)
      println(f,"# nreps: ", nreps )
      println(f,"# max_iterations: ", max_iterations )
      CSV.write( f, df, append=true, writeheader=true )
      println("csvfile written")
    end
  end
  df
end

# Generate phenotypes, sort them by fitness, and return the position of ph in this sorted list of phenotypes
function fitness_rank( p::Parameters, fitness::Vector{Float64}, ph::Goal )::Int64
  fitfunct( g::Goal ) = fitness[ g[1]+1 ]::Float64   # fitness of a Goal phenotype
  fitfunct( g::MyInt ) = fitness[ g+1 ]::Float64   # fitness of a MyInt phenotype
  fitfunct( ch::Circuit ) = fitfunct( output_values(ch) )::Float64  # fitness of a genotype
  #println("ph: ",ph)
  phenotypes = collect(map( MyInt, 0:2^2^p.numinputs-1 ))   # phenotypes are MyInts rather than  Vectors of MyInts 
  sort!( phenotypes, by=x->fitness[x+1], rev=true )
  #println("phenotypes: ",phenotypes)
  #println("fitnesses:  ",map(x->fitfunct(x),phenotypes))
  findfirst(x->x==ph[1],phenotypes)
end
  
# Returns the result of a hill_climb starting from a random chromosome    
# The result is a triple consisting of:
#  1)  The resulting Goal (phenotypes)
#  2)  The fitness of the resulting phenotype
#  3)  The fitness rank of the resulting phenotype where rank 1 is the highest fitness.
function hill_climb_rank( fitness::Vector{Float64}, p::Parameters, funcs::Vector{Func}, max_iterations::Int64)::Tuple{Goal, Float64, Int64 }
  ch = random_chromosome(p,funcs)
  ( hch, hfit, hph) =  hill_climb( fitness, ch, funcs, max_iterations)
  #println("hph: ",hph)
  fitrank = fitness_rank( hch.params, fitness, hph )
  ( hph, hfit, fitrank )
end

# Does nreps repetitions of calls to hill_climb_rank
function run_hill_climb_rank( fitness::Vector{Float64}, p::Parameters, funcs::Vector{Func}, max_iterations::Int64, nreps::Int64; csvfile::String="" )::DataFrame
  df = DataFrame( iteration=Int64[], result_ph=Vector{MyInt}[], fitvalue=Float64[], fitrank=Int64[] )
  fit_goal_list = Tuple{Float64,Goal}[]
  for i = 1:nreps
    ( hph, hfit, fitrank ) = hill_climb_rank( fitness, p, funcs, max_iterations)
    #( hch, hfit, hph) =  hill_climb( fitness, p, funcs, max_iterations)
    row = (i, hph, hfit, fitrank )
    push!(df,row)
  end
  if length(csvfile) > 0
    println("csvfile: ",csvfile)
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs )
      print_parameters(f,p,comment=true)
      println(f,"# nreps: ", nreps )
      println(f,"# max_iterations: ", max_iterations )
      CSV.write( f, df, append=true, writeheader=true )
      println("csvfile written")
    end
  end
  df
end
