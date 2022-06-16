# Test the shape-space property

# Shape-space-covering A GP map has the shape space covering property if, given a phenotype, only a
#  small radius around a sequence encoding that phenotype needs to be explored in order to find the most
#  common phenotypes [Schuster, p, Fontana, Stadler.
#  "From sequences to shapes and back: a case study in RNA secondary structures. Proc R Soc Lond B 1994;255:279–84"].

# This version returns the dataframe of results with one row for each phenotype in goal_list.
# For each phenotype in goal_list, a starting circuit is evolved using function pheno_evolve() which is defined in Evolve.jl
# If pheno_evolve() fails, it returns (nothing,nothing), and no circuit is added to circuits list.
# Thus, there is no row in the output dataframe for this phenotype of goal_list.
function shape_space_counts( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, goal_list::GoalList, max_tries::Int64, max_steps::Int64; 
      use_lincircuit::Bool=false, csvfile::String="" )
  circuits_list = use_lincircuit ? LinCircuit[] : Chromosome[]
  for goal in goal_list
    (c,steps) = pheno_evolve( p, funcs, goal, max_tries, max_steps, use_lincircuit=use_lincircuit )
    if (c,steps) != (nothing,nothing
      push!( circuits_list, c )
    end
  end
  shape_space_counts( p, funcs, num_mutates, circuits_list, use_lincircuit=use_lincircuit, csvfile=csvfile )
end

# This version returns the dataframe of results with one row for each of num_circuits random circuits.
function shape_space_counts( p::Parameters, funcs::Vector{Func}, num_mutates::Int64, num_circuits::Int64; use_lincircuit::Bool=false, csvfile::String="" )
  circuits_list = use_lincircuit ? [ rand_lcircuit( p, funcs ) for _ = 1:num_circuits ] : [ random_chromosome( p, funcs ) for _ = 1:num_circuits ] 
  shape_space_counts( p, funcs, num_mutates, circuits_list, use_lincircuit=use_lincircuit, csvfile=csvfile )
end    

# produces a dataframe with the lengths of the result of pheno_set_funct() for each of num_circuits random circuits.
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
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
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

function genotype_evolvability( c::Circuit, funcs::Vector{Func} )
  outputs_list = mutate_all( c, funcs, output_outputs=true, output_circuits=false )
  length(unique(outputs_list))
end

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

