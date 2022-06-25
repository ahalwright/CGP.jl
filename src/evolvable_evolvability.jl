# Compute phenotype evolvability for a phenotype ph in ph_list by evolving ncircuits circuits that map to ph,
# applying mutate_all() to these circuits, and recording the phenotypes in an count vector indexed # over all phenotypes.
# Produce a dataframe with phenotypes in ph_list as rows and phenotypes mapped to by mutational neighbors of the circuits.
function evolvable_pheno_df( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64; 
    use_lincircuit::Bool=false, csvfile::String="" )
  nrepeats = nprocs()==1 ? 1 : nprocs()-1
  ncircuits_iter = Int(ceil(ncircuits/nrepeats))   # ncircuits/nrepeats rounded up to the next integer.
  result_dict = ph_vect_dict = Dict{Goal,Vector{Int64}}()
  dicts = pmap( i->evolvable_pheno_dict( p, funcs, ph_list, ncircuits_iter, max_tries, max_steps ), 1:nrepeats )
  #dicts = map( i->evolvable_pheno_dictt( p, funcs, ph_list, ncircuits_iter, max_tries, max_steps ), 1:nrepeats )
  result_dict = ph_vect_dict = Dict{Goal,Vector{Int64}}()
  result_dict = merge(+,result_dict,dicts...)
  key_list = sort([k for k in keys(result_dict)])
  df = DataFrame()
  df.pheno_list =  map( k->@sprintf("0x%04x",k[1]), key_list)
  df.evolvability = map( k->length(findall(x->x.!=0,result_dict[k])), key_list )
  df.pheno_vects = map( k->result_dict[k], key_list )
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# ncircuits: ",ncircuits)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function evolvable_pheno_dict( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64; use_lincircuit::Bool=false )
  ph_vect_dict = Dict{Goal,Vector{Int64}}()
  for ph in ph_list
    pheno_count_vect = zeros(Int64,2^2^p.numinputs)   # zero vector over all phenotypes
    evolvable_pheno_count!( pheno_count_vect,p,funcs,ph,ncircuits,max_tries,max_steps) 
    ph_vect_dict[ph] = pheno_count_vect
  end
  ph_vect_dict
end

# Compute phenotype evolvability for a phenotype ph by evolving ncircuits circuits that map to ph,  applying mutate_all() 
#  to these circuits, and recording the phenotypes in a count vector indexed over all phenotypes.
function evolvable_pheno_count!( pheno_count_vect::Vector{Int64}, p::Parameters, funcs::Vector{Func}, ph::Goal, ncircuits::Int64, 
    max_tries::Int64, max_steps::Int64; use_lincircuit::Bool=false )
  for i = 1:ncircuits
    (circ,steps) = pheno_evolve( p, funcs, ph::Goal, max_tries::Int64, max_steps::Int64; use_lincircuit=use_lincircuit )
    (outputs_list,circ_list) = mutate_all( circ, funcs, output_outputs=true, output_circuits=true )
    for pheno in outputs_list
      #println("pheno: ",pheno)
      pheno_count_vect[pheno[1]+1] += 1
    end
  end
end

# Returns a boolean (true/false) matrix of the cases where evolvability succeeds
function pheno_vects_to_boolean_matrix( pheno_vects::Vector{Vector{Int64}} )
  result_matrix = zeros(Bool,length(pheno_vects[1]),length(pheno_vects[1]))
  for i = 1:length(pheno_vects)
    for j = 1:length(pheno_vects[i])
      result_matrix[i,j] = pheno_vects[i][j] != 0 ? true : false
    end
  end
  result_matrix
end

# Calls the previous version when pheno_vects is a vector strings
function pheno_vects_to_boolean_matrix( pheno_vects::Vector{String} )
  ph_vects =map(i->eval(Meta.parse(pheno_vects[i])), 1:length(pheno_vects))  # convert strings to Int64s
  pheno_vects_to_boolean_matrix( ph_vects )
end

# returns a matrix of counts of the phenotypes that contribute to the evolvability set
function pheno_vects_to_evolvable_matrix( pheno_vects::Vector{Vector{Int64}} )
  result_matrix = zeros(Int64,length(pheno_vects[1]),length(pheno_vects[1]))
  for i = 1:length(pheno_vects)
    for j = 1:length(pheno_vects[i])
      result_matrix[i,j] = pheno_vects[i][j] 
    end
  end
  result_matrix
end

# Calls the previous version when pheno_vects is a vector strings
function pheno_vects_to_evolvable_matrix( pheno_vects::Vector{String} )
  ph_vects =map(i->eval(Meta.parse(pheno_vects[i])), 1:length(pheno_vects))  # convert strings to Int64s 
  pheno_vects_to_evolvable_matrix( ph_vects )
end

#  E = pheno_vects_to_evolvable_matrix( pdf.pheno_vects )
function submatrix_to_dataframe( p::Parameters, funcs::Vector{Func}, E::Matrix{Int64}, evdf::DataFrame )
  common_list_int = [0x0000, 0x0005, 0x0011, 0x0022, 0x0033, 0x0044, 0x0055, 0x005f, 0x0077, 0x0088, 0x00a0, 0x00aa, 0x00bb, 0x00cc, 0x00dd, 0x00ee, 0x00fa, 0x00ff]
  common_list_str = ["0x0000","0x0005","0x0011","0x0022","0x0033","0x0044","0x0055","0x005f","0x0077","0x0088","0x00a0","0x00aa","0x00bb","0x00cc","0x00dd","0x00ee","0x00fa","0x00ff"]
  rare_list_str = ["0x0049", "0x0061", "0x0069", "0x006d", "0x0079", "0x0086", "0x0092", "0x0096", "0x009e", "0x00b6"]
  rare_list_int = [ 0x0049, 0x0061, 0x0069, 0x006d, 0x0079, 0x0086, 0x0092, 0x0096, 0x009e, 0x00b6 ]
  common_indices = findall(x->(x in common_list_str),evdf.pheno_list)
  rare_indices = findall(x->(x in rare_list_str),evdf.pheno_list)
  edf = DataFrame()
  for i = 1:length(common_list_str)
    insertcols!(edf,i,Pair(@sprintf("0x%02x",common_indices[i]-1),E[rare_indices,common_indices[i]]))
  end 
  insertcols!(edf,1,"pheno"=>pdf.pheno_list[rare_indices])
  edf
end
