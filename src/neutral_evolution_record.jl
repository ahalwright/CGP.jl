# Evolves a Chromsome or LinCircuit that maps to g starting with chromosome c.
# max_steps is the maximum number of evolutionary steps.
# If evolution hasn't succeeeded in max_steps, return nothing.
# Similar to mut_evolve except that this takes a single goal instead of a goal list as an argument.
# Returns a DataFrame with statistics about the run.
# If onlyShowEpochChange is false, returns information about every neutral or distance-reduncing step
# If onlyShowEpochChange is true, returns information about distance-reduncing step but not neutral steps
#   In this case, returns information about the previous epoch, which is the preceding neutral steps after the last distance-reducing step.
# Evolvability and robustness statistics are genotype evolvability and robustness
function neutral_evolution_record(
  c::Circuit,
  funcs::Vector{Func},
  g::Goal,
  max_steps::Integer;
  print_steps::Bool = false,
  csvfile::String = "",
  onlyShowEpochChange::Bool = false,
) 

  p = c.params
  default_funcs(p)
  #Creating a dictionary to map phenotypes to their complexities
  kdict = kolmogorov_complexity_dict(p, funcs)
  rdict = redundancy_dict(p, funcs)

  println("goal: ",g,"  decimal goal: ",g[1], "   goal Kcomplexity: ", kdict[g[1]] )
  println("max_steps: ",max_steps)
  save_kcomplexities = true  # Save kcomplexity, evolvability, robustness of phenos produced by mutate_all on every step
  if (onlyShowEpochChange)
    save_kcomplexities = false
  end

  # Inner function inside the scope of function neutral_evolution_record()
  # Creates a function within the scope of neutral_evolution and changes the variables without needing them as parameters
  # Can create the inner function anywhere, as long as the function is created prior to its call
  function save_kcomplexity_epochalsteps()
    # ov is the output value of the last genotype of the previous epoch, new_ov is the output value of the new epoch.
    println("function save_kcomplexity_epochalsteps()  step: ", step )
    if step==131 || step==132
      print_circuit( priorEpochGenos_list[end] )
    end

    #Getting new kcomplexity of the new epoch
    push!(newEpochalKComplexity_list, kdict[new_ov[1]])

    #Adds the step count of the epochal step forward
    push!(epochalStep_list, step)
    #Adds the neutral step count of the epochal step forward
    push!(neutralStep_list, neutral_step)

    push!(newPheno_list,new_ov[1])

    #Adds the distance to the goal from the new phenotype
    push!(distanceFromGoal_list, new_distance)

    allPhenosList = MyInt[]  # A list of all phenotypes produced by mutate_all() applied to all genotypes of the prior epoch
    genoPhenosList = Vector{MyInt}[]  # each element is a tuple of the output values and the list of phenos produced by applying mutate_all() to one genotype

    numberTargetFound = reduce( +, map( x -> ((x == g[1]) ? 1 : 0), allPhenosList ))
    push!(numOfTargetFoundInPrevEpoch_list, numberTargetFound)

    # The loop to compute evolvability, lgredund, and robustness  
    mean_robust2 = 0.0
    mean_lg_redund = 0.0   # Establish scope
    if length(genoPhenosList) > 0
      mean_evo_count = reduce( +, map(phl -> (length(phl) > 0 ? length(unique(phl)) - 1 : 0 ), genoPhenosList ))/length(genoPhenosList)
      robust_sum = 0.0
      lg_redund_sum = 0.0
      for i = 1:length(priorEpochGenos_list)
        ovg = output_values(priorEpochGenos_list[i])
        phl = genoPhenosList[i]
        lg_redund = lg10(rdict[ovg[1]])
        lg_redund_sum += lg_redund
        robust = length( filter(z -> ( z==ovg[1] ), phl ) )/length(phl)
        robust_sum += robust
      end
      mean_lg_redund = lg_redund_sum/length(genoPhenosList) 
      mean_robust = robust_sum/length(genoPhenosList)
      mean_robust2 = mean( map( circ->robustness(circ,funcs), priorEpochGenos_list ) )
    else
      mean_lg_redund = 0.0
      mean_evo_count = 0
      mean_robust = 0.0
    end
    println("mean_robust: ",mean_robust, "  mean_robust2: ",mean_robust2)
    push!(averageEvolInPrevEpoch_list, mean_evo_count)
    push!(averageLgRedundInPrevEpoch_list, mean_lg_redund )
    push!(averageRobustInPrevEpoch_list, mean_robust)

    #= Correctness check for evolvability
    genotypeEvolvabilityList = []
    for geno in priorEpochGenos_list
      push!(genotypeEvolvabilityList, length(unique(map(x -> x[1], mutate_all(geno, funcs, output_outputs = true)))))
    end
    mean_evo = length(genotypeEvolvabilityList) > 0 ? mean(genotypeEvolvabilityList) : 0
    println("mean_evo_count: ",mean_evo_count,"  mean_evo: ",mean_evo)
    #push!(averageEvolInPrevEpoch_list, mean_evo)
    =#

    #Calculates the evolvability of the new genotype
    sizeTupleNewGenoEvolv = length(unique(mutate_all(new_c, funcs, output_outputs = true)))
    push!(firstGenoEvol_list, sizeTupleNewGenoEvolv)

    #Calculates the robustness of the new genotype
    println("new_ov: ",new_ov,"  output_values(new_c): ",output_values(new_c))
    phenoListNewGenoRobust = mutate_all(new_c, funcs, output_outputs = true)
    samePhenoCounterNewGenoRobust = 0
    for pheno in phenoListNewGenoRobust
      if (pheno[1] == new_ov[1])
        samePhenoCounterNewGenoRobust += 1
      end
    end
    #println("first robust (samePhenoCounterNewGenoRobust/length(phenoListNewGenoRobust)): ",samePhenoCounterNewGenoRobust/length(phenoListNewGenoRobust),"  robust: ",robustness(new_c,funcs))
    push!( firstGenoRobust_list, (samePhenoCounterNewGenoRobust / length(phenoListNewGenoRobust)))
  end
  # End of inner function save_kcomplexity_epochalsteps()

  LinCirc = typeof(c) == LinCircuit ? :true : :false
  outputPheno_list = MyInt[] #Used to show the current phenotype on every iteration
  step_list = Int64[] # Only used if save_kcomplexities==true
  status_list = String[]  # Only used if save_kcomplexities==true
  kcomplexity_list = Float64[]  # Only used if save_kcomplexities==true
  evolvability_list = Int64[]  # Only used if save_kcomplexities==true
  robust_list = Float64[]   # Only used if save_kcomplexities==true
  lgredund_list = Float64[]   # Only used if save_kcomplexities==true
  intersect_list = Int64[]  # Only used if save_kcomplexities==true 

  #Data to be returned in the case of onlyShowEpochChange being true-------------------------------------------
  #Shows the new phenotype following an ephocal change
  newPheno_list = MyInt[]
  #Shows the step the ephocal change was made on
  epochalStep_list = Int64[]
  neutralStep_list = Int64[]
  #Shows the distance from the new phenotype to the target phenotype
  distanceFromGoal_list = Float64[]
  # Shows K complexity and log redundancy of the new phenotype
  newEpochalKComplexity_list = Float64[]
  newEpochalLGredund_list = Float64[]

  #Shows the amount of times the target phenotype was found in the preceding epoch after a call to mutate all
  numOfTargetFoundInPrevEpoch_list = Int64[]
  #Shows the average redundancy of each genotype in prior epoch
  averageLgRedundInPrevEpoch_list = Float64[]
  #Shows the average evolvability of each genotype in prior epoch
  averageEvolInPrevEpoch_list = Float64[]
  #Records average robustness of each genotype in prior epoch
  averageRobustInPrevEpoch_list = Float64[]

  #The following vectors will record kcomp, lgredund, evolvability and robustness of the first
  #genotype in the new epoch, will also only be used in onlyShowEpochChange == true
  firstGenoEvol_list = Float64[]
  firstGenoRobust_list = Float64[]

  #Vector will store all the genotypes of the "prior" epoch (prior from the POV of a new epoch being disovered)
  #Will be emptied with each new epoch and will be filled with each neutral evolution
  priorEpochGenos_list = Circuit[]

  ov = output_values(c)
  new_c =c  # Establish scope
  new_ov = ov  # Establish scope
  #println("ov: ",ov,"  g: ",g)
  current_distance = hamming_distance(ov[1], g[1], c.params.numinputs)
  new_distance = 0.0  # Value not used.  Establish scope.
  step = 0
  neutral_step = 0   # counts the number of neutral and improvement steps
  while step < max_steps && ov != g   # The main evolutionary loop
    step += 1
    #println("A step: ",step," current_distance: ",current_distance)
    new_c = deepcopy(c)
    if typeof(c) == Chromosome
      (new_c, active) = mutate_chromosome!(new_c, funcs)
    elseif typeof(c) == LinCircuit
      #println("c: ",c)
      new_c = mutate_circuit!(new_c, funcs)
    end
    new_ov = output_values(new_c)
    new_distance = hamming_distance(new_ov, g, c.params.numinputs)
    if step==132
      println( "step: ", step, "  new_ov: ", new_ov, "  new_distance: ", new_distance )
    end
    if new_ov == ov  #Neutral 
      c = new_c
      neutral_step += 1

      #pushes the neutral genotype onto the list
      if onlyShowEpochChange == true
        push!(priorEpochGenos_list, c)
      end

      if save_kcomplexities
        save_kcomplexity(
          new_c,
          c,
          new_ov[1],
          g,
          step,
          "neutral",
          step_list,
          status_list,
          kcomplexity_list,
          evolvability_list,
          robust_list,
          lgredund_list,
          intersect_list,
          outputPheno_list,
          funcs,
          kdict,
        )
      end
      if print_steps && step < 500
        print(
          "step: ",
          step,
          " is pheno neutral.  new_ov: ",
          new_ov,
          "  new_distance: ",
          new_distance,
          "  ",
        )
        print_circuit(new_c)
      end
    elseif new_distance == current_distance #Neutral
      c = new_c
      neutral_step += 1

      #pushes the neutral genotype onto the list
      if onlyShowEpochChange == true
        push!(priorEpochGenos_list, c)
      end

      if save_kcomplexities
        save_kcomplexity(
          new_c,
          c,
          new_ov[1],
          g,
          step,
          "neutral",
          step_list,
          status_list,
          kcomplexity_list,
          evolvability_list,
          robust_list,
          lgredund_list,
          intersect_list,
          outputPheno_list,
          funcs,
          kdict,
        )
      end
      if print_steps && step < 500
        print(
          "step: ",
          step,
          " is fitness neutral.  new_ov: ",
          new_ov,
          "  new_distance: ",
          new_distance,
          "  ",
        )
        print_circuit(new_c)
      end
    elseif new_distance < current_distance   # improvement
      neutral_step += 1
      if print_steps
        println(
          "step: ",
          step,
          "  new_output: ",
          new_ov,
          " distance improved from ",
          current_distance,
          " to ",
          new_distance,
        )
        print_circuit(new_c)
      end
      if save_kcomplexities
        save_kcomplexity(
          new_c,
          c,
          new_ov[1],
          g,
          step,
          "improve",
          step_list,
          status_list,
          kcomplexity_list,
          evolvability_list,
          robust_list,
          lgredund_list,
          intersect_list,
          outputPheno_list,
          funcs,
          kdict,
        )
      end

      #Need to add if onlyShowEpochChange is true to call function
      if (onlyShowEpochChange == true)
        save_kcomplexity_epochalsteps()

        #Empties the list
        empty!(priorEpochGenos_list)

        #adds the new circuit to the priorEpochlist
        push!(priorEpochGenos_list, new_c)
      end
      c = new_c
      ov = new_ov
      current_distance = new_distance
      #println("B step: ",step," current_distance: ",current_distance)
    else
      if print_steps && step <= 20
        print(
          "step: ",
          step,
          "  new_output: ",
          new_ov,
          " current distance: ",
          current_distance,
          " new: ",
          new_distance,
          "  ",
        )
        print_circuit(new_c)
      end
    end
  end # 
  if step == max_steps
    println("neutral evolution failed with ", step, " steps for goal: ", g)
    #return (c, step)
  else
    println("neutral evolution succeeded at step ", step, " for goal: ", g)
    @assert output_values(c) == g
    #return (c, step)
  end
  if save_kcomplexities
    df = DataFrame(
      :step => step_list,
      :status => status_list,
      :kcomplexity => kcomplexity_list,
      :lgredund => lgredund_list,
      :evolvability => evolvability_list,
      :robustness => robust_list,
      :count_goals => intersect_list,
      :pheno => outputPheno_list,
    )
  elseif onlyShowEpochChange
    df = DataFrame(
      :newpheno => newPheno_list,
      :epochalStep => epochalStep_list,
      :neutralStep => neutralStep_list,
      :dist_goal => distanceFromGoal_list,
      :first_geno_evol => firstGenoEvol_list,
      :first_geno_rbst => firstGenoRobust_list,
      :epochal_Kcmplx => newEpochalKComplexity_list,
      #:epochal_lgredund => newEpochalLGredund_list,
      :nTargetFound => numOfTargetFoundInPrevEpoch_list,
      :epochal_evolv => averageEvolInPrevEpoch_list,
      :epochal_lgredund => averageLgRedundInPrevEpoch_list,
      :epochal_rbst => averageRobustInPrevEpoch_list,
    )
  end
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open(csvfile, "w") do f
      println(f, "# date and time: ", Dates.now())
      println(f, "# host: ", hostname, " with ", nprocs() - 1, "  processes: ")
      print_parameters(f, p, comment = true)
      println(f, "# funcs: ", funcs)
      println(f, "# goal: ", g)
      println(f, "# goal K complexity: ", kdict[g[1]])
      println(f, "# max_steps: ", max_steps)
      CSV.write(f, df, append = true, writeheader = true)
    end
  end
  df
end

function save_kcomplexity(
  new_c::Circuit,
  prev_c::Circuit,
  ov::MyInt,
  g::Goal,
  step::Int64,
  status::String,
  step_list::Vector{Int64},
  status_list::Vector{String},
  kcomplexity_list::Vector{Float64},
  evolvability_list::Vector{Int64},
  robust_list::Vector{Float64},
  lgredund_list::Vector{Float64},
  intersect_list::Vector{Int64},
  outputPheno_list::Vector{MyInt},
  funcs::Vector{Func},
  kdict::Dict{MyInt,Int64},
)
  #Adds step count from current run to given list
  push!(step_list, step)

  #Adds the status from the current run to a given list (If the evolution was neutral or improve
  push!(status_list, status)

  #Returns all phenotypes that  result from the given circuit
  phenos = map(x -> x[1], mutate_all(new_c, funcs, output_outputs = true)) # converts goals returned by mutate_all() to MyInts
  #println("phenos: ",phenos)
  #println("ph int g: ",length(findall( x->x==g[1], phenos )))

  #Finds the complexity of the current steps' output values (phenotype)
  push!(kcomplexity_list, kdict[ov])

  #Adds the current mapped phenotype to the list of phenotypes
  push!(outputPheno_list, ov)

  push!(evolvability_list, length(unique(phenos)) - 1)   # Do not include the robustness pheno
  push!(robust_list, length(findall(x -> x == ov, phenos)) / length(phenos))
  push!(lgredund_list, lg10(rdict[ov]))
  push!(intersect_list, length(findall(x -> x == g[1], phenos)))
end
