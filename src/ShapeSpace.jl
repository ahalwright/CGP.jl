# Test the shape-space property
# Shape-space-covering A GP map has the shape space covering property if, given a phenotype, only a
#  small radius around a sequence encoding that phenotype needs to be explored in order to find the most
#  common phenotypes [Schuster, p, Fontana, Stadler.
#  "From sequences to shapes and back: a case study in RNA secondary structures. Proc R Soc Lond B 1994;255:279–84"].

# Runs shape_space_multiple_genos() and returns the fraction of the phenos in phenos_to_cover that are in the output phenotypes of shape_space_multiple_genos().
# If output_phenos==true, then for each 
function fraction_coverage( p::Parameters, funcs::Vector{Func}, phenos_to_cover::Vector{MyInt},  num_mutates::Int64, goal_list::GoalList, circuits_per_goal_list::Vector{Int64},  
      max_tries::Int64, max_steps::Int64; use_lincircuit::Bool=false, output_phenos::Bool=false, csvfile::String="" )
  df = shape_space_multiple_genos( p, funcs, num_mutates, goal_list, circuits_per_goal_list, max_tries, max_steps, output_phenos=true )
  fraction_coverage( df, phenos_to_cover )
end

# df is the dataframe produced by shape_space_multiple_genos() with output_phenos==true
function fraction_coverage( df::DataFrame, phenos_to_cover::Vector{MyInt}; csvfile::String="" )
  n = length(phenos_to_cover)
  fract_covered = zeros(Float64,size(df)[1])
  for i = 1:size(df)[1]
    fract_covered[i] = length( intersect( phenos_to_cover, df.phenos[i] ))/n
  end
  cdf = insertcols!(deepcopy(df), size(df)[2]-1, :fract_covered=>fract_covered )
  cdf
end

# phenotypes that ShapeSpace is trying to cover.  I. e., the "common" phenotypes in the statement of the Shape Space covering hypothesis
# pdf = read_dataframe("../data/7_8_22/evolvable_evolvability_3x1_7_4ch_scmplxP.csv") 
function phenos_to_cover( pdf::DataFrame, quant::Float64 )
  map(x->eval(Meta.parse(x)),pdf[ pdf.ints7_4 .>= quantile(pdf.ints7_4,quant),1])  # Assumes that the first column is the string representation of goals
end
  
# For each phenotype in goal_list, uses function pheno_evolve_to_goals() to evolve a circuit that map to the phenotype.
# For each phenotype, computes the phenotypes in the num_mutates-neighborhoods of this circuit 
# For each phenotype, there is a row of the output dataframe which includes size of these neighborhoods, 
#   and if output_phenos==true, the list of phenos.[
#### Each entry of circuits_per_goal_list is a number of circuits used to compute the unique phenotype count
#function shape_space_multiple_genos( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, goal_list::GoalList, circuits_per_goal_list::Vector{Int64}, 
#      max_tries::Int64, max_steps::Int64; increase_mutates::Bool=false, use_lincircuit::Bool=false, output_phenos::Bool=false, csvfile::String="" )
# in the evolvability paper, circuits_per_goal == 1.
function shape_space_multiple_genos( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, goal_list::GoalList, circuits_per_goal::Int64, 
      max_tries::Int64, max_steps::Int64; use_lincircuit::Bool=false, output_phenos::Bool=false, csvfile::String="" )
  #println("output_phenos: ",output_phenos)
  #println("circuits_per_goal_list: ",circuits_per_goal_list,"  circuits_per_goal_list[end]: ",circuits_per_goal_list[end])
  kdict = kolmogorov_complexity_dict( p, funcs )
  #if circuits_per_goal_list[end] > max_tries
  #  error("circuits_per_goal_list[end] < max_tries.  Increase max_tries in the call to shape_space_multiple_genos() ")
  #end
  #list_of_circuit_steps_lists = pheno_evolve_to_goals( p, funcs, goal_list, circuits_per_goal_list[end], max_tries, max_steps )
  list_of_circuit_steps_lists = pheno_evolve_to_goals( p, funcs, goal_list, circuits_per_goal, max_tries, max_steps )
  #println("typeof(list_of_circuit_steps_lists): ",typeof(list_of_circuit_steps_lists),"  length(list_of_circuit_steps_lists): ",length(list_of_circuit_steps_lists))
  df = DataFrame()
  df.goal=String[] 
  df.complexity=Float64[] 
  df.Kcomplexity=Float64[] 
  df.evolvability=Float64[] 
  df.pheno_counts = Float64[]
  if output_phenos
    df.phenos = Goal[]
  end
  #println("df: ",df,"  names(df): ",names(df))
  # Each call to process_genotype() will generate one row of the dataframe
  df_row_list = pmap( circuit_steps_list->process_genotype( p, funcs, num_mutates, circuits_per_goal, circuit_steps_list, kdict,
      use_lincircuit=use_lincircuit, output_phenos=output_phenos ), list_of_circuit_steps_lists )
  #df_row_list = map( circuit_steps_list->process_genotype( p, funcs, num_mutates, circuits_per_goal, circuit_steps_list, kdict,
  #   use_lincircuit=use_lincircuit, output_phenos=output_phenos ), list_of_circuit_steps_lists )
  #println("size(df): ",size(df))
  for df_row in df_row_list
    #println("length(df_row): ",length(df_row))
    push!(df,df_row)
  end
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      #println(f,"# circuits_per_goal_list: ",circuits_per_goal_list)
      println(f,"# circuits_per_goal: ",circuits_per_goal)
      println(f,"# num_mutates: ",num_mutates)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      #println(f,"# length(circuits_list): ",length(circuits_list))
      CSV.write( f, df, append=true, writeheader=true )
    end
    df
  end
  df
end

# helper function for function shape_space_multiple_genos()
function process_genotype( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, circuits_per_goal::Int64,
      circuit_steps_list::Union{Vector{Tuple{LinCircuit,Int64}},Vector{Tuple{Chromosome,Int64}}}, kdict::Dict{MyInt,Int64}; 
    use_lincircuit::Bool=false, output_phenos::Bool=false )
  #println("output_phenos: ",output_phenos)
  #kdict = kolmogorov_complexity_dict( p, funcs )
  #circuits_list = map( cs->cs[1], circuits_steps_list )  # doesn't work for some reason
  circuits_list = [ circuit_steps_list[i][1] for i = 1:length(circuit_steps_list) ]
  #steps_list = map( cs->cs[2], circuits_steps_list )
  steps_list = [ circuit_steps_list[i][2] for i = 1:length(circuit_steps_list) ] 
  # convert from circuits to circuit_ints for pheno_set_funct()
  circuit_int_list = map(c->(use_lincircuit ? circuit_to_circuit_int(c,funcs) : chromosome_to_int(c)), circuits_list )
  pheno = output_values(circuits_list[1])
  #println("pheno: ",pheno)
  pheno_string = @sprintf("0x%04x",pheno[1])
  #println("pheno_string: ",pheno_string)
  complexity_list = map( complexity5, circuits_list )
  complexity = sum(complexity_list)/length(complexity_list)
  Kcomplexity_list = map( ph->kdict[ph], map(c->output_values(c)[1],circuits_list ))
  Kcomplexity = sum(Kcomplexity_list)/length(Kcomplexity_list)
  evolvability_list = map( c->genotype_evolvability(c,funcs), circuits_list )
  evolvability = sum(evolvability_list)/length(evolvability_list)
  #current_pheno(ci,p,funcs) = output_values( use_lincircuit ? circuit_int_to_circuit(ci,p,funcs) : int_to_chromosome(ci,p,funcs) )
  cumm_pheno_set = Set(Goal[])
  #pheno_counts = zeros(Int64,circuits_per_goal)
  pheno_counts = 0
  println("[circuit_int_list[1:circuits_per_goal]][1]: ",[circuit_int_list[1:circuits_per_goal]][1])
  new_pheno_set = pheno_set_funct( p, funcs, [circuit_int_list[1:circuits_per_goal]][1], num_mutates, Set([pheno]) )
  union!( cumm_pheno_set, new_pheno_set )
  pheno_counts += length( cumm_pheno_set )
  #println("pheno_counts: ",pheno_counts)
  if output_phenos
    df_row = vcat([ pheno_string, complexity, Kcomplexity, evolvability ], pheno_counts )
    cumm_pheno_list = [ ph[1] for ph in cumm_pheno_set ]   # Assumes 1 output by extracting 1 component
    push!(df_row,cumm_pheno_list)
  else
    df_row = [ pheno_string, complexity, Kcomplexity, evolvability, pheno_counts ]
  end
  #println("df_row: ",df_row)
  return df_row 
end

# This version returns the dataframe of results with one row for each phenotype in goal_list.
# For each phenotype in goal_list, a starting circuit is evolved using function pheno_evolve() which is defined in Evolve.jl
# If pheno_evolve() fails, it returns (nothing,nothing), and no circuit is added to circuits list.
# Thus, there is no row in the output dataframe for this phenotype of goal_list.
function shape_space_counts( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, goal_list::GoalList, max_tries::Int64, max_steps::Int64; 
      use_lincircuit::Bool=false, csvfile::String="" )
  circuits_list = use_lincircuit ? LinCircuit[] : Chromosome[]
  for goal in goal_list
    (c,steps) = pheno_evolve( p, funcs, goal, max_tries, max_steps, use_lincircuit=use_lincircuit )
    if (c,steps) != (nothing,nothing)
      push!( circuits_list, c )
    end
  end
  #println("length(circuits_list): ",length(circuits_list))
  shape_space_counts( p, funcs, num_mutates, circuits_list, use_lincircuit=use_lincircuit, csvfile=csvfile )
end

# This version returns the dataframe of results with one row for each of num_circuits random circuits.
function shape_space_counts( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, num_circuits::Int64; use_lincircuit::Bool=false, csvfile::String="" )
  #println("num_circuits: ",num_circuits)
  circuits_list = use_lincircuit ? [ rand_lcircuit( p, funcs ) for _ = 1:num_circuits ] : [ random_chromosome( p, funcs ) for _ = 1:num_circuits ] 
  shape_space_counts( p, funcs, num_mutates, circuits_list, use_lincircuit=use_lincircuit, csvfile=csvfile )
end    

# Write a DataFrame with one row per phenotype.  
# DataFrame has 3 columns:  :goal, :pheno_counts, and :pheno.  
# :pheno column is the phenotypes discovered by doing num_mutates mutations of a circuit that maps to the corresponding phenotype/.
# circuit_ints_list_list is the circuits_list column from a count output file such as those in data/counts.
# goals_list is a list of MyInt phenotypes that will be processed.
# Example:  wdf = read_dataframe("../data/counts/count_outputs_ch_4funcs_3inputs_8gates_4lb_W.csv")
# Note that p must be the same parameters as those for the counts file.
# shape_space_circuit_ints_list( p, funcs, 1, wdf.circuits_list )
function shape_space_circuit_ints_list( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, circuit_ints_list_list::Union{Vector{Vector{Int128}},Vector{String}}; 
     goals_list::Vector{MyInt}=MyInt[], use_lincircuit::Bool=false, csvfile::String="" )
  if length(goals_list) > 0
    goals_indices = map(x->x+MyInt(1),goals_list)   # convert to 1-based indexing
    circuit_ints_list_list = circuit_ints_list_list[goals_indices]   # restrict to goals_indices
  end
  if typeof(circuit_ints_list_list) == Vector{String}
    # Convert elements of circuits_ints_list_list from strings to Vector{Int128}
    circuit_ints_list_list = map( x->eval(Meta.parse(x)), circuit_ints_list_list )
  end
  goal_list = Vector{MyInt}[]
  pheno_counts_list = Int64[]
  phenos_list_list = Vector{MyInt}[]
  results = pmap( circuit_ints_list->shape_space_circuit_ints( p, funcs, num_mutates, circuit_ints_list), circuit_ints_list_list )
  #results = map( circuit_ints_list->shape_space_circuit_ints( p, funcs, num_mutates, circuit_ints_list), circuit_ints_list_list )
  for res in results
    (goal,phenos) = res
    push!( goal_list, goal )
    #println("goal: ",goal,"  goal_list: ",goal_list)
    push!(pheno_counts_list, length(phenos))
    push!( phenos_list_list, map(x->x[1], phenos) )
  end
  if map( x->x[1], goal_list) != goals_list
    println("goal_list not the same as goals_list in function shape_space_circuit_ints_list() ")
    println("goal_list: ",map( x->x[1], goal_list))
    println("goals_list: ",goals_list)
  end
  df = DataFrame( :goal=>goal_list, :pheno_count=>pheno_counts_list, :phenos=>phenos_list_list )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# num_mutates: ",num_mutates)
      println(f,"# length(circuit_ints_list_list): ",length(circuit_ints_list_list))
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Produce a goal phenotype and list of phenos from an element of the circuit_ints field of a counts file
# Helper function for function shape_space_circuit_ints_list()
function shape_space_circuit_ints( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, circuit_ints_list::Vector{Int128};
  use_lincircuit::Bool=false )
  default_funcs(p)   # Set global variable Ones to be appropriate for parameters p.
  # Define a local function that computes the phenotype from a circuit_int ci.
  current_pheno(ci,p,funcs) = output_values( use_lincircuit ? circuit_int_to_circuit(ci,p,funcs) : int_to_chromosome(ci,p,funcs) )
  circuits_list = use_lincircuit ? map( ci->int_to_circuit( ci, p, funcs ), circuit_ints_list) : map( ci->int_to_chromosome( ci, p, funcs ), circuit_ints_list)
  outputs_list = map( output_values, circuits_list )
  #println("outputs_list[1]: ",outputs_list[1])
  if length( unique(outputs_list) )  != 1 
    error("outputs_list includes more than 1 phenotype in function shape_space_circuit_ints().  Check that parameters agrees with dataframe for circuit_ints_list")
  end
  goal = outputs_list[1]
  ci = rand( circuit_ints_list )
  phenos = set_to_list(pheno_set_funct( p, funcs, [ci], num_mutates, Set([current_pheno(ci,p,funcs)]) ))
  (goal,phenos)
end  

# produces a dataframe with the lengths of the result of pheno_set_funct() for each the circuits in circuits_list.
function shape_space_counts( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, circuits_list::Union{ Vector{LinCircuit}, Vector{Chromosome} }; 
    use_lincircuit::Bool=false, csvfile::String="" )
  #circuits_list = use_lincircuit ? [ rand_lcircuit( p, funcs ) for _ = 1:num_circuits ] : [ random_chromosome( p, funcs ) for _ = 1:num_circuits ] 
  circuit_int_list = unique(map(c->(use_lincircuit ? circuit_to_circuit_int(c,funcs) : chromosome_to_int(c)), circuits_list ))
  phenos_list = map( c->output_values(c), circuits_list )
  println("phenos_list: ",phenos_list)
  pheno_strings_list = map( ph->@sprintf("0x%04x",ph[1]), phenos_list )
  #println("pheno_strings_list: ",pheno_strings_list)
  complexity_list = map( complexity5, circuits_list )
  evolvability_list = map( c->genotype_evolvability(c,funcs), circuits_list )
  current_pheno(ci,p,funcs) = output_values( use_lincircuit ? circuit_int_to_circuit(ci,p,funcs) : int_to_chromosome(ci,p,funcs) )
  #result_list = map( ci->length(pheno_set_funct( p, funcs, [ci], num_mutates, Set([current_pheno(ci,p,funcs)]) )), circuit_int_list )
  result_list = pmap( ci->length(pheno_set_funct( p, funcs, [ci], num_mutates, Set([current_pheno(ci,p,funcs)]) )), circuit_int_list )
  df = DataFrame()
  df.phenos = pheno_strings_list
  df.complexity = complexity_list
  df.evolvability = evolvability_list
  df.results = result_list
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# num_mutates: ",num_mutates)
      println(f,"# length(circuits_list): ",length(circuits_list))
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

    
function pheno_evolve_to_goals( p::Parameters, funcs::Vector{Func}, goal_list::GoalList, num_circuits_per_goal::Int64, max_tries::Int64, max_steps::Int64;
      use_lincircuit::Bool=false )
  list_of_circuits_steps_lists = pmap( goal->pheno_evolve( p, funcs, goal, num_circuits_per_goal, max_tries, max_steps, use_lincircuit=use_lincircuit ), goal_list)
  #list_of_circuits_steps_lists = map( goal->pheno_evolve( p, funcs, goal, num_circuits_per_goal, max_tries, max_steps, use_lincircuit=use_lincircuit ), goal_list)
end

function genotype_evolvability( c::Circuit, funcs::Vector{Func} )
  outputs_list = mutate_all( c, funcs, output_outputs=true, output_circuits=false )
  length(unique(outputs_list))
end

# Recursively compute the set of phenotypes reached by all num_mutates mutations from all circuits in circuit_int_list
function pheno_set_funct( p::Parameters, funcs::Vector{Func}, circuit_int_list::Vector{Int128}, num_mutates::Int64, pheno_set::Set{Vector{MyInt}};
    use_lincircuit::Bool=false )::Set{Vector{MyInt}}
  println("num_mutates: ",num_mutates,"  length(circuit_int_list): ",length(circuit_int_list),"  length(pheno_set): ",length(pheno_set))
  #println("pheno_set: ",sort([x[1] for x in pheno_set]))
  circuit_list = use_lincircuit ? map( ci->circuit_int_to_circuit(ci,p,funcs), circuit_int_list) : map( ci->int_to_chromosome(ci,p,funcs), circuit_int_list )
  if num_mutates > 0
    new_circ_ints_list = Int128[]
    for ci in circuit_int_list
      circuit = use_lincircuit ? circuit_int_to_circuit(ci,p,funcs) : int_to_chromosome(ci,p,funcs)
      (outputs_list,circ_list) = mutate_all( circuit, funcs, output_outputs=true, output_circuits=true )
      circ_ints_list = unique(map(c->(use_lincircuit ? circuit_to_circuit_int(c,funcs) : chromosome_to_int(c)), circ_list ))
      #println("outputs_list: ",outputs_list)
      union!( pheno_set, Set(outputs_list) )  
      union!(new_circ_ints_list,circ_ints_list) 
    end
    return pheno_set_funct( p, funcs, new_circ_ints_list, num_mutates-1, pheno_set )
  else
    return pheno_set
  end
end

# Returns a list of the robustnesses of numcircuits random genotypes
function genotype_robustness( p::Parameters, funcs::Vector{Func}, pheno::Goal, numcircuits::Int64; use_lincircuit::Bool=false )
  robust_list = Float64[]
  for i = 1:numcircuits
    c = use_lincircuit ? rand_lcircuit( p, funcs ) : random_chromosome( p, funcs )
    push!( robust_list, robustness( c, funcs ) )
  end
  robust_list
end

# list of phenotypes in ph_list  which are above the quantile_value quantile for ph_list
function common_phenos( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, quantile_value::Float64, rdict::Dict{MyInt,Int64} )
  #rdict = redundancy_dict(p,funcs)
  redundancy_list = map( ph->lg10(rdict[ph[1]]), ph_list )
  #println("redundancy_list: ",redundancy_list)
  cutoff_value = quantile( redundancy_list, quantile_value )
  println("cutoff_value: ",cutoff_value)
  common_list = filter( ph->lg10(rdict[ph[1]]) >= cutoff_value, ph_list )
end
  
# computes the fraction of common phenotypes covered by num_mutates mutations of the given phenotype.
# ph_df should be the DataFrame result of function shape_space_multiple_genos() or function shape_space_circuit_ints_list() with option output_phenos=true.
# Note that num_mutates is a parameter of shape_space_multiple_genos() and of function shape_space_circuit_ints_list()
# There is no check that parameters p and funcs agree with the settings used to generate DataFrame ph_df.
function shape_space_fract_successes( p::Parameters, funcs::Vector{Func}, quantile_values::Union{Float64,Vector{Float64}}, ph_df::DataFrame; csvfile::String="" )
  rdict = redundancy_dict(p,funcs)
  if !( "phenos" in names(ph_df) )
    error("ph_df must include a column \"phenos\" in function shape_space_fraction")
  end
  if typeof(ph_df.phenos[1]) == String
    ph_df.phenos = map( x->eval(Meta.parse(x)), ph_df.phenos )
  end
  if typeof(ph_df.goal[1]) <: AbstractString
    ph_df.goal = map( x->[eval(Meta.parse(x))], ph_df.goal )
  #else
    #ph_df.goal = map( x->[x], ph_df.goal )
  end
  println("ph_df.goal: ",ph_df.goal)
  if typeof(quantile_values) <: AbstractFloat
    quantile_values = [quantile_values]   # convert to list
  end
  all_phenos_list = map(x->[x],MyInt(0):MyInt(2^2^p.numinputs-1))  # all phenotypes for these parameters
  phenoSet_list = map( i->Set(ph_df.phenos[i]),1:size(ph_df)[1])
  for qv in quantile_values
    commonPhenos = common_phenos( p, funcs, all_phenos_list, qv, rdict )
    numcommon = length(commonPhenos)
    println("quantile_values: ",qv,"  numcommon: ",numcommon)
    commonPhenosSet = Set( map(x->x[1], commonPhenos ) )
    setdiff_list = map( phset->setdiff(commonPhenosSet,phset), phenoSet_list )
    insert_position = size(ph_df)[2]  
    insertcols!(ph_df,insert_position,Symbol("fs$(qv)")=>map(x->(numcommon-length(x))/numcommon,setdiff_list))
  end
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# quantile_values: ",quantile_values)
      println(f,"# size(ph_df): ",size(ph_df))
      #println(f,"# max_tries: ",max_tries)
      #println(f,"# max_steps: ",max_steps)
      CSV.write( f, ph_df, append=true, writeheader=true )
    end
  end
  ph_df
end

# circuit_ints_list_list is the circuits_list column from a count output file such as those in data/counts.
# goals_list is a list of MyInt phenotypes that will be processed.
# Example:  wdf = read_dataframe("../data/counts/count_outputs_ch_4funcs_3inputs_8gates_4lb_W.csv")
# Note that p must be the same parameters as those for the counts file.
# shape_space_sampling_successes( p, funcs, 2, [0.9,0.95], wdf.circuits_list )
function shape_space_sampling_successes( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, quantile_values::Vector{Float64}, circuit_ints_list_list::Union{Vector{Vector{Int128}},Vector{String}};
     goals_list::Vector{MyInt}=MyInt[], use_lincircuit::Bool=false, csvfile::String="" )
  ph_df = shape_space_circuit_ints_list( p, funcs, num_mutates, circuit_ints_list_list, goals_list=goals_list, use_lincircuit=use_lincircuit )
  println("ph_df.goal: ",ph_df.goal)
  ph_df = shape_space_fract_successes( p, funcs, quantile_values, ph_df )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nworkers(),"  processes: " )
      println(f,"# sampling")
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# num_mutates: ", num_mutates)
      println(f,"# quantile_values: ",quantile_values)
      println(f,"# size(ph_df): ",size(ph_df))
      CSV.write( f, ph_df, append=true, writeheader=true )
    end
  end
  ph_df
end

function shape_space_evolution_successes( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, quantile_values::Vector{Float64}, ph_int_list::Vector{MyInt}, max_tries::Int64, max_steps::Int64;
     use_lincircuit::Bool=false, circuits_per_goal::Int64=1, csvfile::String="" )
  goal_list = map( ph->[ph], ph_int_list )
  ph_df = shape_space_multiple_genos( p, funcs, num_mutates, goal_list, circuits_per_goal, max_tries, max_steps, use_lincircuit=use_lincircuit, output_phenos=true )
  ph_df = shape_space_fract_successes( p, funcs, quantile_values, ph_df )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nworkers(),"  processes: " )
      println(f,"# evolution")
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs)
      println(f,"# max_tries: ", max_tries)
      println(f,"# max_steps: ", max_steps)
      println(f,"# num_mutates: ", num_mutates)
      println(f,"# quantile_values: ",quantile_values)
      println(f,"# size(ph_df): ",size(ph_df))
      CSV.write( f, ph_df, append=true, writeheader=true )
    end
  end
  ph_df
end
