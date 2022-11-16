# Test the shape-space property
# Shape-space-covering A GP map has the shape space covering property if, given a phenotype, only a
#  small radius around a sequence encoding that phenotype needs to be explored in order to find the most
#  common phenotypes [Schuster, p, Fontana, Stadler.
#  "From sequences to shapes and back: a case study in RNA secondary structures. Proc R Soc Lond B 1994;255:279–84"].

# Runs shape_space_multiple_genos() and returns the fraction of the phenos in phenos_to_cover that are in the output phenotypes of shape_space_multiple_genos().
# If output_phenos==true, then for each 
function fraction_coverage( p::Parameters, funcs::Vector{Func}, phenos_to_cover::Vector{MyInt},  num_mutates::Int64, goal_list::GoalList, circuits_per_goal_list::Vector{Int64},  
      max_tries::Int64, max_steps::Int64; increase_mutates::Bool=false, use_lincircuit::Bool=false, output_phenos::Bool=false, csvfile::String="" )
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
  
# For each phenotype in goal_list, uses function pheno_evolve_to_goals() to evolve circuits that map to the phenotype.
# For each phenotype, computes the phenotypes in the unions of the k-neighborhoods of these circuits for each k in circuits_per_goal_list[2:end].  k==num_mutates.
# For each phenotype, there is a row of the output dataframe which includes size of these neighborhoods for each k, 
#   and if output_phenos==true, the list of phenos for the maximum k value.
# Each entry of circuits_per_goal_list is a number of circuits used to compute the unique phenotype count
function shape_space_multiple_genos( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, goal_list::GoalList, circuits_per_goal_list::Vector{Int64}, 
      max_tries::Int64, max_steps::Int64; increase_mutates::Bool=false, use_lincircuit::Bool=false, output_phenos::Bool=false, csvfile::String="" )
  println("circuits_per_goal_list: ",circuits_per_goal_list,"  circuits_per_goal_list[end]: ",circuits_per_goal_list[end])
  kdict = kolmogorov_complexity_dict( p, funcs )
  if circuits_per_goal_list[end] > max_tries
    error("circuits_per_goal_list[end] < max_tries.  Increase max_tries in the call to shape_space_multiple_genos() ")
  end
  list_of_circuit_steps_lists = pheno_evolve_to_goals( p, funcs, goal_list, circuits_per_goal_list[end], max_tries, max_steps )
  println("typeof(list_of_circuit_steps_lists): ",typeof(list_of_circuit_steps_lists),"  length(list_of_circuit_steps_lists): ",length(list_of_circuit_steps_lists))
  df = DataFrame()
  df.goal=String[] 
  df.complexity=Float64[] 
  df.Kcomplexity=Float64[] 
  df.evolvability=Float64[] 
  #println("names(df): ",names(df),"  size(df): ",size(df))
  num_pheno_cols = increase_mutates ? num_mutates : length(circuits_per_goal_list)-1
  pheno_cols = [ Symbol("pheno_count$(i)") for i = 1: num_pheno_cols ]
  for phc in pheno_cols
    df[:,phc] = Float64[]
  end
  if output_phenos
    df.phenos = Goal[]
  end
  #println("df: ",df,"  names(df): ",names(df))
  # Each call to process_genotype() will generate one row of the dataframe
  #df_row_list = pmap( circuit_steps_list->process_genotype( p, funcs, num_mutates, goal_list, circuits_per_goal_list, circuit_steps_list, 
  df_row_list = pmap( circuit_steps_list->process_genotype( p, funcs, num_mutates, goal_list, circuits_per_goal_list, circuit_steps_list, 
      increase_mutates=increase_mutates, use_lincircuit=use_lincircuit, output_phenos=output_phenos ), list_of_circuit_steps_lists )
  #df_row_list = map( circuit_steps_list->process_genotype( p, funcs, num_mutates, goal_list, circuits_per_goal_list, circuit_steps_list, 
  #    increase_mutates=increase_mutates, use_lincircuit=use_lincircuit, output_phenos=output_phenos ), list_of_circuit_steps_lists )
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
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# circuits_per_goal_list: ",circuits_per_goal_list)
      println(f,"# increase_mutates: ",increase_mutates)
      println(f,"# num_mutates: ",num_mutates)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      #println(f,"# length(circuits_list): ",length(circuits_list))
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function process_genotype( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, goal_list::GoalList, circuits_per_goal_list::Vector{Int64},
      circuit_steps_list::Union{Vector{Tuple{LinCircuit,Int64}},Vector{Tuple{Chromosome,Int64}}}; 
    increase_mutates::Bool=false, use_lincircuit::Bool=false, output_phenos::Bool=false )
  kdict = kolmogorov_complexity_dict( p, funcs )
  #circuits_list = map( cs->cs[1], circuits_steps_list )  # doesn't work for some reason
  circuits_list = [ circuit_steps_list[i][1] for i = 1:length(circuit_steps_list) ]
  #steps_list = map( cs->cs[2], circuits_steps_list )
  steps_list = [ circuit_steps_list[i][2] for i = 1:length(circuit_steps_list) ] 
  circuit_int_list = map(c->(use_lincircuit ? circuit_to_circuit_int(c,funcs) : chromosome_to_int(c)), circuits_list )
  pheno = output_values(circuits_list[1])
  println("pheno: ",pheno)
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
  pheno_counts = zeros(Int64,length(circuits_per_goal_list)-1)
  if !increase_mutates
    for j = 2:length(circuits_per_goal_list)
      #println("j: ",j,"  [circuit_int_list[j-1:j]]: ",[circuit_int_list[circuits_per_goal_list[j-1]:circuits_per_goal_list[j]]][1])
      new_pheno_set = pheno_set_funct( p, funcs, [circuit_int_list[circuits_per_goal_list[j-1]:circuits_per_goal_list[j]]][1],num_mutates, Set([pheno]) )
      union!( cumm_pheno_set, new_pheno_set )
      pheno_counts[j-1] = length( cumm_pheno_set )
    end
  else
    for num_mutate = 1:num_mutates
      #println("j: ",j,"  [circuit_int_list[j-1:j]]: ",[circuit_int_list[circuits_per_goal_list[j-1]:circuits_per_goal_list[j]]][1])
      new_pheno_set = pheno_set_funct( p, funcs, [circuit_int_list[end]][1], num_mutate, Set([pheno]) )  # Doesn't work 6/22/22
      union!( cumm_pheno_set, new_pheno_set )
      pheno_counts[num_mutate-1] = length( cumm_pheno_set )
    end
  end
  #println("pheno_counts: ",pheno_counts)
  if output_phenos
    df_row = vcat([ pheno_string, complexity, Kcomplexity, evolvability ], pheno_counts )
    cumm_pheno_list = [ ph[1] for ph in cumm_pheno_set ]   # Assumes 1 output by extracting 1 component
    push!(df_row,cumm_pheno_list)
  else
    df_row = vcat([ pheno_string, complexity, Kcomplexity, evolvability ], pheno_counts )
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
  #result_list = map( ci->(pheno_set_funct( p, funcs, [ci], num_mutates, Set([current_pheno(ci,p,funcs)]) )), circuit_int_list )
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
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# num_mutates: ",num_mutates)
      println(f,"# length(circuits_list): ",length(circuits_list))
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end
    
#=
function pheno_evolve_to_goals( p::Parameters, funcs::Vector{Func}, goal_list::GoalList, num_circuits_per_goal::Int64, max_tries::Int64, max_steps::Int64;
      use_lincircuit::Bool=false )
  list_of_circuits_steps_lists = map( goal->pheno_evolve( p, funcs, goal, num_circuits_per_goal, max_tries, max_steps, use_lincircuit=use_lincircuit ), GoalList )
end
=#

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
function pheno_set_funct( p::Parameters, funcs::Vector{Func}, circuit_int_list::Vector{Int128}, num_mutates::Int64, pheno_set::Set;
    use_lincircuit::Bool=false )
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

