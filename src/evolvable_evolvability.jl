# Compute phenotype evolvability for a phenotype ph in ph_list by either evolving ncircuits circuits that map to ph,
#  or by converting circ_int_list to circuits.  Tnen mutate_all() is applied to these circuits, and a count vector of phenotype counts is returned..
# Produce a dataframe with phenotypes in ph_list as rows and phenotypes mapped to by mutational neighbors of the circuits.
function evolvable_pheno_df( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64; circ_int_lists::Vector{Vector{Int128}}=Vector{Int128}[], use_lincircuit::Bool=false, csvfile::String="" )
  nrepeats = nprocs()==1 ? 1 : nprocs()-1
  ncircuits_iter = Int(ceil(ncircuits/nrepeats))   # ncircuits/nrepeats rounded up to the next integer.
  dict_pairs = pmap( i->evolvable_pheno_dict( p, funcs, ph_list, ncircuits_iter, max_tries, max_steps, circ_int_lists=circ_int_lists, use_lincircuit=use_lincircuit ), 1:nrepeats )
  #dict_pairs = map( i->evolvable_pheno_dict( p, funcs, ph_list, ncircuits_iter, max_tries, max_steps, circ_int_lists=circ_int_lists, use_lincircuit=use_lincircuit ), 1:nrepeats )
  ph_dicts = map(x->x[1],dict_pairs)
  cmplx_dicts = map(x->x[2],dict_pairs)
  result_ph_dict = ph_vect_dict = Dict{Goal,Vector{Int64}}()
  result_ph_dict = merge(+,result_ph_dict,ph_dicts...)   # Merge the ph_dicts 
  result_cmplx_dict = Dict{Goal,Vector{Float64}}()
  result_cmplx_dict = merge(vcat,result_cmplx_dict,cmplx_dicts...) # Merge the cmplx_dicts 
  #println("result_cmplx_dict: ",result_cmplx_dict)
  ph_key_list = sort([k for k in keys(result_ph_dict)])
  cmplx_key_list = sort([k for k in keys(result_cmplx_dict)])
  df = DataFrame()
  df.pheno_list =  map( k->@sprintf("0x%04x",k[1]), ph_key_list)
  df.evolvability = map( k->length(findall(x->x.!=0,result_ph_dict[k])), ph_key_list )
  df.complexity = map( k->sum(result_cmplx_dict[k])/length(result_cmplx_dict[k]), cmplx_key_list )
  df.pheno_vects = map( k->result_ph_dict[k], ph_key_list )
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# ncircuits: ",ncircuits)
      println(f,"# length(circ_int_lists): ",length(circ_int_lists))
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function evolvable_pheno_dict( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64; 
    circ_int_lists::Vector{Vector{Int128}}, use_lincircuit::Bool=false )
  ph_vect_dict = Dict{Goal,Vector{Int64}}()
  ph_cmplx_dict = Dict{Goal,Vector{Float64}}()
  for i = 1:length( ph_list )
    pheno_count_vect = zeros(Int64,2^2^p.numinputs)   # zero vector over all phenotypes
    ph = ph_list[i]
    if length(circ_int_lists) > 0
      circuit_list = compute_circuit_list( p, funcs, ph, ncircuits, max_tries, max_steps, circ_int_list=circ_int_lists[i], use_lincircuit=use_lincircuit )
      evolvable_pheno_count!( pheno_count_vect, circuit_list, funcs, circ_int_list=circ_int_lists[i], use_lincircuit=use_lincircuit ) 
    else
      circuit_list = compute_circuit_list( p, funcs, ph, ncircuits, max_tries, max_steps, circ_int_list=Int128[], use_lincircuit=use_lincircuit )
      evolvable_pheno_count!( pheno_count_vect, circuit_list, funcs, circ_int_list=Int128[], use_lincircuit=use_lincircuit ) 
    end
    complexity_list = use_lincircuit ? map(c->lincomplexity(c,funcs), circuit_list)  : map(c->complexity5(c), circuit_list)
    #println("complexity_list: ",complexity_list)
    ph_cmplx_dict[ph] = complexity_list
    ph_vect_dict[ph] = pheno_count_vect
  end
  (ph_vect_dict, ph_cmplx_dict)
end

function compute_circuit_list( p::Parameters, funcs::Vector{Func}, ph::Goal, ncircuits::Int64,
    max_tries::Int64, max_steps::Int64; circ_int_list::Vector{Int128}, use_lincircuit::Bool=false )
  if length(circ_int_list) == 0
    circuit_list = use_lincircuit ? LinCircuit[] : Chromosome[]
    for i = 1:ncircuits
      (circ,steps) = pheno_evolve( p, funcs, ph, max_tries, max_steps; use_lincircuit=use_lincircuit )
      push!(circuit_list,circ)
    end
  else
    circuit_list = map( ci->( use_lincircuit ? (circuit_int_to_circuit( ci, p, funcs )) : (int_to_chromosome( ci, p, funcs ))), circ_int_list )
  end 
  circuit_list
end

# Compute phenotype evolvability for a phenotype ph by either evolving ncircuits circuits that map to ph or by using the given circ_int_list,
#  then applying mutate_all() to each circuit in circ_int_list.
#  The vector pheno_count_vect, which is indexed over phenotypes, is modified in place
function evolvable_pheno_count!( pheno_count_vect::Vector{Int64}, circuit_list::Union{Vector{Chromosome},Vector{LinCircuit}}, funcs::Vector{Func}; 
    circ_int_list::Vector{Int128}, use_lincircuit::Bool=false )
  for circ in circuit_list
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

# For CGP 3x1 7 gates 4 lb:
 CGP_common_list_int = [0x0000, 0x0005, 0x0011, 0x0022, 0x0033, 0x0044, 0x0055, 0x005f, 0x0077, 0x0088, 0x00a0, 0x00aa, 0x00bb, 0x00cc, 0x00dd, 0x00ee, 0x00fa, 0x00ff]
 CGP_common_list_str = ["0x0000","0x0005","0x0011","0x0022","0x0033","0x0044","0x0055","0x005f","0x0077","0x0088","0x00a0","0x00aa","0x00bb","0x00cc","0x00dd","0x00ee","0x00fa","0x00ff"]
 CGP_rare_list_str = ["0x0049", "0x0061", "0x0069", "0x006d", "0x0079", "0x0086", "0x0092", "0x0096", "0x009e", "0x00b6"]
 CGP_rare_list_int = [ 0x0049, 0x0061, 0x0069, 0x006d, 0x0079, 0x0086, 0x0092, 0x0096, 0x009e, 0x00b6 ]
 LGP_rare_list_str=["0x0016", "0x0061", "0x0068", "0x0069", "0x006b", "0x006d", "0x0086", "0x0092", "0x0094", "0x0096", "0x00b6", "0x00d6", "0x00e9"]  # Complexities >= 3.6
 LGP_common_list_str=["0x0000","0x0001","0x0002","0x0004","0x0010","0x0011","0x003f","0x007f","0x0080","0x00df","0x00ef","0x00f7","0x00fb","0x00fc","0x00fe","0x00ff"] # complexities <= 2.6
#  E = pheno_vects_to_evolvable_matrix( pdf.pheno_vects )
#  Possible problem if the evdf dataframe is directly generated rather than read from a CSV file.
# Note that default values are from common_str and rare_str.
function submatrix_to_dataframe( p::Parameters, funcs::Vector{Func}, E::Matrix{Int64}, evdf::DataFrame, common_list::Union{Vector{MyInt},Vector{String}}=common_str, 
    rare_list::Union{Vector{MyInt},Vector{String}}=rare_str; source::String="rare", dest::String="common" )
  println("source: ",source,"  dest: ",dest)
  to_string(x::MyInt) = MyInt==UInt16 ? @sprintf("0x%02x",x) : @sprintf("0x%04x",x)
  common_list_str = typeof(common_list)==Vector{String} ? common_list : map(x->to_string(x), common_list )
  rare_list_str = typeof(rare_list)==Vector{String} ? rare_list : map(x->to_string(x), rare_list )
  common_indices = findall(x->(x in common_list_str),evdf.pheno_list)
  rare_indices = findall(x->(x in rare_list_str),evdf.pheno_list)
  println("rare_indices: ",rare_indices)
  println("common_indices: ",common_indices)
  edf = DataFrame()
  for i = 1:(dest=="rare" ? length(rare_list_str) : length(common_list_str))
    #println("i: ",i)
    if source=="rare"
      if dest=="rare"
        insertcols!(edf,i,Pair(@sprintf("0x%02x",rare_indices[i]-1),E[rare_indices,rare_indices[i]]))
      else
        insertcols!(edf,i,Pair(@sprintf("0x%02x",common_indices[i]-1),E[rare_indices,common_indices[i]]))
      end
    else
      if dest=="rare"
        insertcols!(edf,i,Pair(@sprintf("0x%02x",rare_indices[i]-1),E[common_indices,rare_indices[i]]))
      else
        insertcols!(edf,i,Pair(@sprintf("0x%02x",common_indices[i]-1),E[common_indices,common_indices[i]]))
      end
    end
  end 
  if source=="rare"
    insertcols!(edf,1,"pheno"=>evdf.pheno_list[rare_indices])
  else
    insertcols!(edf,1,"pheno"=>evdf.pheno_list[common_indices])
  end
  edf
end

function total_evolvability( pdf::DataFrame )
  to_binary(x::Bool) = x ? 1 : 0
  to_bool(x::Int64) = x != 0 ? true : false
  E = pheno_vects_to_evolvable_matrix( pdf.pheno_vects )
  n = length(pdf.pheno_list)
  total_evo = Int64[]
  for k = 1:n
    bool_vect = map( v->to_bool(v), E[k,:]+E[:,k] )
    evo_vect = map( v->to_binary(v), bool_vect )
    evo_vect[k] = 0
    push!(total_evo,sum(evo_vect))
  end
  total_evo
end

function total_evol( pdf::DataFrame )
  to_binary(x::Bool) = x ? 1 : 0
  to_bool(x::Int64) = x != 0 ? true : false
  B = map( to_bool, pheno_vects_to_evolvable_matrix( pdf.pheno_vects ))  # Boolean evolvability matrix
  map( x->sum( B[x,:] .|| B[:,x] ) - to_binary( B[x,x] ), 1:length(pdf.pheno_vects ) )
end
