# Eves a Chromsome or LinCircuit that maps to g starting with chromosome c.
# max_steps is the maximum number of evolutionary steps.
# If evolution hasn't succeeeded in max_steps, return nothing.
# insert_gate_prob is the probability of inserting a gate on a mutation of a chromosome.
# delete_gate_prob is similar for deleting a gate.
# Similar to mut_evolve except that this takes a single goal instead of a goal list as an argument.
function neutral_evolution_record( c::Circuit, funcs::Vector{Func}, g::Goal, max_steps::Integer; print_steps::Bool=false, save_kcomplexities::Bool=true, csvfile::String="", onlyShowEpochalChange::Bool=false ) # Save kcomplexity, evolvability, robustness of phenos produced by mutate_all on every step
  insert_gate_prob = 0.0
  delete_gate_prob = 0.0
  #if typeof(c) == Chromosome
  #  println("neutral evolution: ")
  #  print_circuit(c)
  #end
  #default_funcs(c.params)
  p = c.params

  #Creating a dictionary to map phenotypes to their complexities
  kdict = kolmogorov_complexity_dict(p,funcs)

  LinCirc = typeof(c) == LinCircuit ? :true : :false
  outputPheno_list = MyInt[] #Used to show the current phenotype on every iteration
  step_list = Int64[] # Only used if save_kcomplexities==true
  status_list = String[]  # Only used if save_kcomplexities==true
  kcomplexity_list = Float64[]  # Only used if save_kcomplexities==true
  evolvability_list = Int64[]  # Only used if save_kcomplexities==true
  robust_list = Float64[]   # Only used if save_kcomplexities==true
  intersect_list = Int64[]  # Only used if save_kcomplexities==true 
  #println("LinCirc: ",LinCirc,"  Ones: ",@sprintf("0x%08x",Ones),"  CGP.Ones: ",@sprintf("0x%08x",CGP.Ones))
  #println("numgates: ",c.params.numinteriors)
  #
  
  #Data to be returned in the case of onlyShowEpochalChange being true
  #Shows the new phenotype following an ephocal change
  newPheno_list = MyInt[]
  #Shows the step the ephocal change was made on
  epochalStep_list = Int64[]
  #Shows the distance from the new phenotype to the target phenotype
  distanceFromGoal_list = Float64[]

  #Shows the amount of times the target phenotype was found in the preceding epoch
  numOfTargetFoundInPrevEpoch_list = Int64[]
  #Shows average kcomplexity of each genotype in prior epoch
  averageKCompInPrevEpoch_list = Float64[]
  #Shows the average redundancy of each genotype in prior epoch
  averageRedundInPrevEpoch_list = Float64[]
  #Shows the average evolvability of each genotype in prior epoch
  averageEvolInPrevEpoch_list = Float64[]
  #Records average robustness of each genotype in prior epoch
  avergaeRobustInPrevEpoch_list = Float64[]

  #The following vectors will record kcomp, redund, evolvability and robustness of the first
  #genotype in the new epoch, will also only be used in onlyShowEpochalChange == true
  firstGenoKComp_list = Float64[]
  firstGenoRedund_list = Float64[]
  firstGenoEvol_list = Float64[]
  firstGenoRobust_list = Float64[]

  #Vector will store all the genotypes of the "prior" epoch (prior from the POV of a new epoch being disovered)
  #Will be emptied with each new epoch and will be filled with each neutral evolution
  priorEpochGenos_list = Circuit[]

  step = 0
  ov = output_values( c )
  #println("ov: ",ov,"  g: ",g)
  current_distance = hamming_distance( ov, g, c.params.numinputs )
  while step < max_steps && ov != g
    step += 1
    #println("A step: ",step," current_distance: ",current_distance)
    new_c = deepcopy(c)
    if typeof(c) == Chromosome
      (new_c,active) = mutate_chromosome!( new_c, funcs, insert_gate_prob=insert_gate_prob, delete_gate_prob=delete_gate_prob )
    elseif typeof(c) == LinCircuit
      #println("c: ",c)
      new_c = mutate_circuit!( new_c, funcs )
    end
    new_ov = output_values( new_c )
    new_distance = hamming_distance( new_ov, g, c.params.numinputs )
    #=
    if step < 20
      print("step: ",step,"  ov: ",ov,"  new_ov: ",new_ov,"  cur dis: ",current_distance,"  new_dis: ",new_distance,"  " )
      print_circuit(new_c)
    end
    =#
    #print("new_c: ")
    #print_circuit(new_c)
    #print("    c: ")
    #print_circuit(c)
    if new_ov == ov  #Neutral 
      c = new_c

      #pushes the neutral genotype onto the list
        #onlyShowEpochalChange
      if onlyShowEpochalChange 
	      push!(priorEpochGenos_list, c)
      end

      if save_kcomplexities
        save_kcomplexity( new_c, c, new_ov[1], g, step, "neutral", step_list, status_list, kcomplexity_list, evolvability_list, robust_list, intersect_list, outputPheno_list, funcs, kdict)
      end
      if print_steps && step < 500
        print("step: ",step," is pheno neutral.  new_ov: ",new_ov,"  new_distance: ",new_distance,"  ")
        print_circuit(new_c)
      end
    elseif new_distance == current_distance #Neutral
      c = new_c

      #pushes the neutral genotype onto the list
      if onlyShowEpochalChange == true
	      push!(priorEpochgenos_list, c)
      end

      if save_kcomplexities
        save_kcomplexity( new_c, c, new_ov[1], g, step, "neutral", step_list, status_list, kcomplexity_list, evolvability_list, robust_list, intersect_list, outputPheno_list, funcs, kdict)
      end
      if print_steps && step < 500
        print("step: ",step," is fitness neutral.  new_ov: ",new_ov,"  new_distance: ",new_distance,"  ")
        print_circuit(new_c)
      end
    elseif new_distance < current_distance   # improvement
      if print_steps
        println("step: ",step,"  new_output: ",new_ov," distance improved from ",current_distance," to ",new_distance)
        print_circuit(new_c)
      end
      if save_kcomplexities
        save_kcomplexity( new_c, c, new_ov[1], g, step, "improve", step_list, status_list, kcomplexity_list, evolvability_list, robust_list, intersect_list, outputPheno_list, funcs, kdict )
      end
      c = new_c
      ov = new_ov
      current_distance = new_distance
      #println("B step: ",step," current_distance: ",current_distance)
    else
      if print_steps && step <= 20
        print("step: ",step,"  new_output: ",new_ov," current distance: ",current_distance," new: ",new_distance,"  ")
        print_circuit(new_c)
      end 
    end
  end
  if save_kcomplexities
    df = DataFrame( :step=>step_list, :status=>status_list, :kcomplexity=>kcomplexity_list, :evolvability=>evolvability_list, :robustness=>robust_list,
        :count_goals=>intersect_list, :output_values=>map(x->[x],outputPheno_list))
  elseif step == max_steps
    println("neutral evolution failed with ",step," steps for goal: ",g)
    return (c, step)
  else
    println("neutral evolution succeeded at step ",step," for goal: ",g)
    @assert output_values(c) == g
    return (c, step)
  end
  if save_kcomplexities  && length(csvfile) > 0
	  hostname = readchomp("hostname")
	  open(csvfile, "w") do f
		  println(f,"# date and time: ",Dates.now())
		  println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
		  print_parameters(f,p,comment=true)
		  println(f,"# funcs: ", funcs)
		  println(f,"# max_steps: ",max_steps)
		  CSV.write(f, df, append=true, writeheader=true )
	  end
  end
  df
end



function save_kcomplexity( new_c::Circuit, prev_c::Circuit, ov::MyInt, g::Goal, step::Int64, status::String, step_list::Vector{Int64},status_list::Vector{String}, kcomplexity_list::Vector{Float64}, evolvability_list::Vector{Int64}, robust_list::Vector{Float64}, intersect_list::Vector{Int64}, pheno_list::Vector{MyInt}, funcs::Vector{Func}, kdict::Dict{MyInt, Int64} )
	  #Adds step count from current run to given list
	  push!( step_list, step )

	  #Adds the status from the current run to a given list (If the evolution was neutral or improve
	  push!( status_list, status )

	  #Returns all phenotypes that  result from the given circuit
	  #The latest code reviewer (11/2/22) Does not know why a map call is needed
	  phenos = map( x->x[1], mutate_all( new_c, funcs, output_outputs=true ) )
	  #println("phenos: ",phenos)
          #println("ph int g: ",length(findall( x->x==g[1], phenos )))
	  
	  #Finds the complexity of the current steps' output values (phenotype)
	  push!( kcomplexity_list, kdict[ov])
	  
	  #Adds the current mapped phenotype to the list of phenotypes
	  push!(pheno_list, ov)

	  #push!( kcomplexity_list, mutual_information( phenos, fill( g[1], length(phenos ) ) ) )
	  push!( evolvability_list, length(unique(phenos) ) )
	  push!( robust_list, length( findall( x->x==ov, phenos ) )/length(phenos) )
          push!( intersect_list, length(findall( x->x==g[1], phenos )))
end



#Function called after epochal success, fills the lists with values to be put into lists
#=
function save_kcomplexity_epochalsteps(new_c::Circuit, current_Pheno::MyInt, kdict::Dict{Myint, Int64}, g::Goal, step::Int64, priorGenos::Vector{Circuit}, newPheno_list::Vector{MyInt}, epochalStep_list::Vector{Int64}, distanceFromGoal_list::Vector{Float64}, numOfTargetFoundInPrevEpoch_list::Vector{Int64}, averageKCompInPrevEpoch_list::Vector{Float64}, averageRedundInPrevEpoch_list::Vector{Float64}, averageEvolInPrevEpoch_list::Vector{Float64}, averageRobustInPrevEpoch_list::Vector{Float64}, firstGenoKComp_list::Vector{Float64}, firstGenoRedund_list::Vector{Float64}, firstGenoEvol_list::Vector{Float64}, firstGenoRobust_list::Vector{Float64})
	#Adds the step count of the epochal step forward
	push!(epochalStep_list, step)

	#Adds the new phenotype to the 

=#



