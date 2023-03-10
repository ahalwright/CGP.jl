# Compute phenotype evolvability for a phenotype ph in ph_list by either evolving ncircuits circuits that map to ph,
#  or by converting circ_int_list to circuits.  Tnen mutate_all() is applied to these circuits, and a count vector of phenotype counts is returned..
# Produce a dataframe with phenotypes in ph_list as rows and phenotypes mapped to by mutational neighbors of the circuits.
using Base.Threads

# epochal evolution version
# if csvfile == "" returns the approximate phenonet adjacency matrix
# if csvfile != "" returns a dataframe representing the approximate phenonet adjacency matrix and saves the dataframe to the
function phenonet_matrix_evol_approx( p::Parameters, funcs::Vector{Func}, ncircuits::Int64, max_steps::Int64, max_tries::Int64; csvfile::String="" )
  phlist = map(x->[x],collect(MyInt(0):MyInt(2^2^p.numinputs-1)))
  pdf = evolvable_pheno_df( p, funcs, phlist, ncircuits, max_tries, max_steps )
  E = pheno_vects_to_evolvable_matrix( pdf.pheno_vects)  # phenonet matrix
  if length(csvfile) > 0
    phl = collect(MyInt(0):MyInt(2^2^p.numinputs-1))
    df = matrix_to_dataframe( E, phl )
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# evolution pheno matrix")
      println(f,"# ncircuits: ",ncircuits)
      println(f,"# length(phl): ",length(phl))
      #println(f,"# length(circ_int_lists): ",length(circ_int_lists))
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# random_walk sampling version
# if csvfile == "" returns the approximate phenonet adjacency matrix
# if csvfile != "" returns a dataframe representing the approximate phenonet adjacency matrix and saves the dataframe to the
function phenonet_matrix_sampling_approx( p::Parameters, funcs::Vector{Func}, nwalks::Int64, steps::Int64; csvfile::String="" )
  phl = collect(MyInt(0):MyInt(2^2^p.numinputs-1))
  rrw_df = run_random_walks_parallel( p, nwalks, phl, steps, exclude_zero_rows=false, output_dict=false, save_complex=false )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# sampling pheno matrix")
      println(f,"# length(phl): ",length(phl))
      println(f,"# nwalks: ",nwalks)
      println(f,"# steps: ",steps)
      CSV.write( f, rrw_df, append=true, writeheader=true )
    end
    rrw_df
  else
    E = df_to_matrix_mt( rrw_df, 5 )  # phenonet matrix
  end
end

# Evolves evolvability, Tononi complexity, and pheno_vects for the phenotypes in ph_list using ncircuits per phenotype
# Evolution version
# Alternative function:  run_geno_complexity() in Evolvability.jl
function evolvable_pheno_df( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64;  use_lincircuit::Bool=false, csvfile::String="" )
  println("evolvable_pheno_df: length(ph_list): ",length(ph_list))
  nrepeats = nprocs()==1 ? 1 : nprocs()-1
  ncircuits_iter = Int(ceil(ncircuits/nrepeats))   # ncircuits/nrepeats rounded up to the next integer.
  dict_pairs = pmap( i->evolvable_pheno_dict( p, funcs, ph_list, ncircuits_iter, max_tries, max_steps, use_lincircuit=use_lincircuit ), 1:nrepeats )
  #dict_pairs = map( i->evolvable_pheno_dict( p, funcs, ph_list, ncircuits_iter, max_tries, max_steps, use_lincircuit=use_lincircuit ), 1:nrepeats )
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
    hostname = readchomp(`hostname`) 
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# ncircuits: ",ncircuits)
      println(f,"# length(ph_list): ",length(ph_list))
      #println(f,"# length(circ_int_lists): ",length(circ_int_lists))
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Computes evolvability, Tononi complexity for the phenotypes in ph_list using ncircuits per phenotype
# Requires circ_int_lists from sampling.
# Sampling version
function evolvable_pheno_df( p::Parameters, funcs::Vector{Func}, phlist::Vector, circ_int_lists::Vector; use_lincircuit::Bool=false, csvfile::String="" )
  if typeof(phlist[1]) <: AbstractString
    phlist = map(x->[eval(Meta.parse(x))], phlist )
  end
  if typeof(circ_int_lists[1]) <: AbstractString
    circ_int_lists = map(x->eval(Meta.parse(x)), circ_int_lists )
  end
  (result_ph_dict,result_cmplx_dict) = evolvable_pheno_dict( p, funcs, phlist, 0, 0, 0, circ_int_lists=circ_int_lists, use_lincircuit=use_lincircuit )
  ph_key_list = sort([k for k in keys(result_ph_dict)])
  cmplx_key_list = sort([k for k in keys(result_cmplx_dict)])
  df = DataFrame()
  df.pheno_list =  map( k->@sprintf("0x%04x",k[1]), ph_key_list)
  df.evolvability = map( k->length(findall(x->x.!=0,result_ph_dict[k])), ph_key_list )
  df.complexity = map( k->sum(result_cmplx_dict[k])/length(result_cmplx_dict[k]), cmplx_key_list )
  df.pheno_vects = map( k->result_ph_dict[k], ph_key_list )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`) 
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      #println(f,"# ncircuits: ",ncircuits)
      #println(f,"# length(ph_list): ",length(ph_list))
      println(f,"# length(circ_int_lists): ",length(circ_int_lists))
      #println(f,"# max_tries: ",max_tries)
      #println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Computes pheno vects for phenotypes in ph_list
# if circ_int_lists is empty, use evolution to return pheno_vects
# if circ_int_lists is nonempty, use it to return pheno_vects
# returns a pair of Dicts:  
#   ph_vect_dict maps phenotypes to pheno_vects
#   ph_cmplx_dict maps phenotypes to Tononi complexities
function evolvable_pheno_dict( p::Parameters, funcs::Vector{Func}, ph_list::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64; 
    circ_int_lists::Vector{Vector{Int128}}=Vector{Int128}[], use_lincircuit::Bool=false )
  default_funcs(p)
  ph_vect_dict = Dict{Goal,Vector{Int64}}()
  ph_cmplx_dict = Dict{Goal,Vector{Float64}}()
  for i = 1:length( ph_list )
    pheno_count_vect = zeros(Int64,2^2^p.numinputs)   # zero vector over all phenotypes
    ph = ph_list[i]
    if length(circ_int_lists) > 0
      circuit_list = compute_circuit_list( p, funcs, ph, 0, 0, 0, circ_int_list=circ_int_lists[i], use_lincircuit=use_lincircuit )
      evolvable_pheno_count!( pheno_count_vect, circuit_list, funcs, circ_int_list=circ_int_lists[i], use_lincircuit=use_lincircuit ) 
    else
      circuit_list = compute_circuit_list( p, funcs, ph, ncircuits, max_tries, max_steps, circ_int_list=Int128[], use_lincircuit=use_lincircuit )
      evolvable_pheno_count!( pheno_count_vect, circuit_list, funcs, circ_int_list=Int128[], use_lincircuit=use_lincircuit ) 
    end
    for circ in circuit_list
      ov = output_values(circ)
      if ov != ph
        println("i: ",i,"  ph: ",ph," ov: ",ov)
      end
    end
    complexity_list = use_lincircuit ? map(c->lincomplexity(c,funcs), circuit_list)  : map(c->complexity5(c), circuit_list)
    #println("complexity_list: ",complexity_list)
    ph_cmplx_dict[ph] = complexity_list
    ph_vect_dict[ph] = pheno_count_vect
  end
  (ph_vect_dict, ph_cmplx_dict)
end

# Computes a circuit list either by:
#   if length(circ_int_list) == 0 uses pheno_evolve() to evolve ncircuits genotypes to map to the given phenotype ph
#   if length(circ_int_list) > 0 converts the circuit ints in circ_int_list to circuits
function compute_circuit_list( p::Parameters, funcs::Vector{Func}, ph::Goal, ncircuits::Int64,
    max_tries::Int64, max_steps::Int64; circ_int_list::Vector{Int128}, use_lincircuit::Bool=false )
  if length(circ_int_list) == 0
    circuit_list = use_lincircuit ? LinCircuit[] : Chromosome[]
    for i = 1:ncircuits
      (circ,steps) = pheno_evolve( p, funcs, ph, max_tries, max_steps; use_lincircuit=use_lincircuit )
      if circ !== nothing
        push!(circuit_list,circ)
      end
    end
  else
    circuit_list = map( ci->( use_lincircuit ? (circuit_int_to_circuit( ci, p, funcs )) : (int_to_chromosome( ci, p, funcs ))), circ_int_list )
  end 
  circuit_list
end

# Compute phenotype evolvability for a phenotype ph by applying mutate_all() to each circuit in circ_int_list.
#  The vector pheno_count_vect, which is indexed over phenotypes, is modified in place
# not multithreaded
function evolvable_pheno_count!( pheno_count_vect::Vector{Int64}, circuit_list::Union{Vector{Chromosome},Vector{LinCircuit}}, funcs::Vector{Func}; 
    circ_int_list::Vector{Int128}, use_lincircuit::Bool=false )
  for circ in circuit_list
    (outputs_list,circ_list) = mutate_all( circ, funcs, output_outputs=true, output_circuits=true )
    for pheno in outputs_list
      # println("pheno: ",pheno,"  Ones: ",Ones)
      pheno_count_vect[pheno[1]+1] += 1
    end
  end
end

# Returns a boolean (true/false) matrix of the cases where evolvability succeeds
# not multithreaded
function pheno_vects_to_boolean_matrix( pheno_vects::Vector{Vector{Int64}} )
  result_matrix = zeros(Bool,length(pheno_vects[1]),length(pheno_vects[1]))
  for i = 1:length(pheno_vects)
    for j = 1:length(pheno_vects[i])
      result_matrix[i,j] = pheno_vects[i][j] != 0 ? 1 : 0
    end
  end
  result_matrix
end

# Calls the previous version when pheno_vects is a vector of strings
function pheno_vects_to_boolean_matrix( pheno_vects::Vector{String} )
  ph_vects =map(i->eval(Meta.parse(pheno_vects[i])), 1:length(pheno_vects))  # convert strings to Int64s
  pheno_vects_to_boolean_matrix( ph_vects )
end

# returns a matrix of counts of the phenotypes that contribute to the evolvability set
# Example from notes/9_5_22.txt:
# pdf = read_dataframe("../data/7_8_22/evolvable_evolvability_3x1_7_4ch_scmplxP.csv")
# E = pheno_vects_to_evolvable_matrix( pdf.pheno_vects )
function pheno_vects_to_evolvable_matrix( pheno_vects::Vector{Vector{Int64}} )
  result_matrix = zeros(Int64,length(pheno_vects[1]),length(pheno_vects[1]))
  for i = 1:length(pheno_vects)
    for j = 1:length(pheno_vects[i])
      result_matrix[i,j] = pheno_vects[i][j] 
    end
  end
  result_matrix
end

# Calls the previous version when pheno_vects is a vector of strings
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
function submatrix_to_dataframe( p::Parameters, funcs::Vector{Func}, E::Matrix{Int64}, evdf::DataFrame, source_list::Union{Vector{MyInt},Vector{String}}, 
    dest_list::Union{Vector{MyInt},Vector{String}}  )
  tostring(x::MyInt) = MyInt==UInt16 ? @sprintf("0x%04x",x) : @sprintf("0x%04x",x)
  source_list_str = typeof(source_list)==Vector{String} ? source_list : map(x->tostring(x), source_list )
  dest_list_str = typeof(dest_list)==Vector{String} ? dest_list : map(x->tostring(x), dest_list )
  #println("dest_list_str: ",dest_list_str)
  source_indices = findall(x->(x in source_list_str),evdf.pheno_list)
  dest_indices = findall(x->(x in dest_list_str),evdf.pheno_list)
  #println("dest_indices: ",dest_indices)
  #println("source_indices: ",source_indices)
  println("sum E submatrix: ",sum(E[source_indices,dest_indices]))
  edf = DataFrame()
  for i = 1:length(dest_list) 
    #println("i: ",i)
    insertcols!(edf,i,@sprintf("0x%02x",source_indices[i]-1)=>E[source_indices,dest_indices[i]]) 
  end 
  insertcols!(edf,1,"pheno"=>evdf.pheno_list[source_indices])
  edf
end

function submatrix_to_dataframe( p::Parameters, funcs::Vector{Func}, E::Matrix{Int64}, evdf::DataFrame, row_col_list::Union{Vector{MyInt},Vector{String}} )
  tostring(x::MyInt) = MyInt==UInt16 ? @sprintf("0x%04x",x) : @sprintf("0x%04x",x)
  row_col_list_str = typeof(row_col_list)==Vector{String} ? row_col_list : map(x->tostring(x[1]), row_col_list )
  bv = BitVector( map( i->evdf.pheno_list[i] in row_col_list_str, 1:size(evdf)[1] ))
  row_col_indices = map(i->findfirst(x->row_col_list_str[i]==x, evdf.pheno_list),1:length(row_col_list_str))
  println("row_col_indices: ",row_col_indices)
  edf = DataFrame()
  for i = 1:length(row_col_list)
    #println("i: ",i)
    insertcols!(edf,i,@sprintf("0x%02x",row_col_indices[i]-1)=>E[row_col_indices,row_col_indices[i]])
  end
  insertcols!(edf,1,"pheno"=>evdf.pheno_list[row_col_indices])
  edf

end

#  println("sum E submatrix: ",sum(E[row_col_indices,dest_indices]))
#  map(i->findfirst(x->row_col_list_str[i]==x, evdf.pheno_list),1:length(row_col_list_str))
function submatrix_count( p::Parameters, funcs::Vector{Func}, E::Matrix{Int64}, evdf::DataFrame, source_list::Union{Vector{MyInt},Vector{String}},
    dest_list::Union{Vector{MyInt},Vector{String}}  )
  tostring(x::MyInt) = MyInt==UInt16 ? @sprintf("0x%04x",x) : @sprintf("0x%04x",x)
  source_list_str = typeof(source_list)==Vector{String} ? source_list : map(x->tostring(x), source_list )
  dest_list_str = typeof(dest_list)==Vector{String} ? dest_list : map(x->tostring(x), dest_list )
  #println("dest_list_str: ",dest_list_str)
  source_indices = findall(x->(x in source_list_str),evdf.pheno_list)
  dest_indices = findall(x->(x in dest_list_str),evdf.pheno_list)
  #println("dest_indices: ",dest_indices)
  #println("source_indices: ",source_indices)
  sum(E[source_indices,dest_indices])
end

# Returns a DataFrame representing a submatrix of the evolvability matrix correcponding to phdf.pheno_vects. 
# The matrix has rows corresponding to source which must be "rare" or "common", and columns corresponding to dest which also must be "rare" or "common"
# The rare and common strings can be defined either on the basis of redundancy or on the basis of k_complexity.  
# In the former case, redund_symbol should be name of the redundancy column of phdf (such as :ints7_4), and in the former case this :k_complexity.
# If 0<common_q<1 and 0<rare_q<1 these are interpreted as redundancy quantiles
# If common_q and rare_q are Int64's, these are interpreted as Kolomogorov complexity cutoffs
# test setup: p = Parameters(3,1,8,5); funcs=default_funcs(p); redund_symbol=:ints8_5; rare_q=0.06; common_q=0.96; source="common"; dest="rare"
# test setup: p = Parameters(3,1,8,5); funcs=default_funcs(p); redund_symbol=:k_complexity; rare_q=6; common_q=1; source="common"; dest="rare"
#   phdf = read_dataframe("../data/7_17_22/evolvable_evolabilityCGP_3x1_8_5_7_17_22cmplxA.csv")
function phenodf_to_submatrix_dataframe( p::Parameters, funcs::Vector{Func}, phdf::DataFrame, redund_symbol::Symbol, rare_q::Number, common_q::Number, source::String, dest::String;
      Bool_output::Bool=false )  # Output matrix is Boolean, i. e, 0/1
  println("source: ",source,"  dest: ",dest)
  if !( source in ["rare","common"]) || !(dest in ["rare","common"])
    error("Error in function phenodf_to_submatrix_dataframe():  both rare and common must be one of the strings 'rare' or 'common'.")
  end
  tostring(x::MyInt) = @sprintf("0x%04x",x)
  #tostring(x::MyInt) = @sprintf("0x%02x",x) 
  E = Bool_output ? pheno_vects_to_boolean_matrix( phdf.pheno_vects ) : pheno_vects_to_evolvable_matrix( phdf.pheno_vects )
  if 0.0 < common_q < 1.0 && 0.0 < rare_q < 1.0   # Use common_q and rare_q as redundancy quantiles
    common = phdf[phdf[:,redund_symbol].>=quantile(phdf[:,redund_symbol],common_q),:pheno_list]
    println("common: ",common)
    rare = phdf[phdf[:,redund_symbol].<=quantile(phdf[:,redund_symbol],rare_q),:pheno_list]
  elseif typeof(common_q) == Int64 && typeof(rare_q) == Int64  # Use common_q and rare_q as Kolmogorov complexity cutoffs
    common = phdf[phdf[:,redund_symbol].<=common_q,:pheno_list]
    println("common: ",common)
    rare = phdf[phdf[:,redund_symbol].>=rare_q,:pheno_list]
  #else
  end
  common_indices = findall(x->(x in common),phdf.pheno_list)
  rare_indices = findall(x->(x in rare),phdf.pheno_list)  
  common_list_str = typeof(common)==Vector{MyInt} ? map(x->tostring(x), common ) : common
  rare_list_str = typeof(rare)==Vector{MyInt} ? map(x->tostring(x), rare ) : rare
  common_indices = findall(x->(x in common_list_str),phdf.pheno_list)
  rare_indices = findall(x->(x in rare_list_str),phdf.pheno_list)
  common_list_str = map(s-> occursin("0x00",s) ? string("0x",s[5:end]) : s, common )  # E.g.: Convert "0x0078" to "0x78"
  rare_list_str = map(s-> occursin("0x00",s) ? string("0x",s[5:end]) : s, rare )  # E.g.: Convert "0x006a" to "0x6a"
  println("rare_list_str: ",rare_list_str)
  println("rare_indices: ",rare_indices)
  println("common_indices: ",common_indices)
  edf = DataFrame( :goals=>source=="rare" ? rare_list_str : common_list_str )
  for i = 1:(dest=="rare" ? length(rare_indices) : length(common_indices))
    #println("i: ",i)
    if source=="rare"
      if dest=="rare"
        insertcols!(edf,i+1,Pair(@sprintf("0x%02x",rare_indices[i]-1),E[rare_indices,rare_indices[i]]))
      else
        insertcols!(edf,i+1,Pair(@sprintf("0x%02x",common_indices[i]-1),E[rare_indices,common_indices[i]]))
      end
    else
      if dest=="rare"
        insertcols!(edf,i+1,Pair(@sprintf("0x%02x",rare_indices[i]-1),E[common_indices,rare_indices[i]]))
      else
        insertcols!(edf,i+1,Pair(@sprintf("0x%02x",common_indices[i]-1),E[common_indices,common_indices[i]]))
      end
    end
  end 
  edf
end

# For the 3x1 7gts 4lb case, pdf = read_dataframe("../data/7_8_22/evolvable_evolvability_3x1_7_4ch_scmplxP.csv")
#  E = pheno_vects_to_evolvable_matrix( pdf.pheno_vects )
#  B is the Boolean verson of matrix E
# Close to the corresponding function in random_walk.jl
function total_evol( pdf::DataFrame )
  to_binary(x::Bool) = x ? 1 : 0
  to_bool(x::Int64) = x != 0 ? true : false
  B = map( to_bool, pheno_vects_to_evolvable_matrix( pdf.pheno_vects ))  # Boolean evolvability matrix
  map( x->sum( B[x,:] .|| B[:,x] ) - to_binary( B[x,x] ), 1:length(pdf.pheno_vects ) )
end

# less elegant version
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

function df_to_matrix( df::DataFrame, startcolumn::Int64; denormalize::Bool=false )  
  matrix = zeros(Float64,size(df)[1],size(df)[2]-startcolumn+1)
  for i = 1:size(df)[1]
    for j = startcolumn:size(df)[2]
      #println("i,j: ",(i,j))
      matrix[i,j-startcolumn+1] = df[i,j]
    end
  end
  if denormalize
    if !("redund" in names(phdf))
      println("dataframe must include a column :redund in function df_to_matrix")
    else
      for i = 1:size(matrix)[1]
        if !isnan(phdf.redund[i])
          matrix[i,:] *= phdf.redund[i]
        end
      end
    end
  end
  matrix
end

# Multi-threaded version
function df_to_matrix_mt( df::DataFrame, startcolumn::Int64; denormalize::Bool=false )  
  matrix = zeros(Float64,size(df)[1],size(df)[2]-startcolumn+1)
  #for i = 1:size(df)[1]
  Threads.@threads for i = 1:size(df)[1]
    for j = startcolumn:size(df)[2]
      #println("i,j: ",(i,j))
      matrix[i,j-startcolumn+1] = df[i,j]
    end
  end
  if denormalize
    if !("redund" in names(phdf))
      println("dataframe must include a column :redund in function df_to_matrix")
    else
      for i = 1:size(matrix)[1]
        if !isnan(phdf.redund[i])
          matrix[i,:] *= phdf.redund[i]
        end
      end
    end
  end
  matrix
end
# Compute multiset average of the entropy evolvability of a list of phenotypes.
# Based on evolving goals to map to the phenotypes

function run_multiset_evolvability( p::Parameters, funcs::Vector{Func}, phlist::GoalList, nreps::Int64, max_tries::Int64, max_steps::Int64;
    csvfile::String="" )
  lg10(x) = x==0 ? 0 : log10(x)
  k_complex_dict = kolmogorov_complexity_dict( p )
  k_complex_list = map(x->k_complex_dict[x[1]],phlist)
  redund_dict = redundancy_dict( p )
  redund_list = map(x->lg10(redund_dict[x[1]]),phlist)
  entropy_list = Float64[]
  len_nonzeros_list = Int64[]
  for ph in phlist
    @sprintf("0x%04x",ph[1])
    (pcounts, entropy, len_nonzeros ) = multiset_phenotype_evolvability( p, funcs, ph, nreps, max_tries, max_steps )
    push!(entropy_list,entropy)
    push!(len_nonzeros_list,len_nonzeros)
  end
  df = DataFrame( :goal=>map(x->@sprintf("0x%04x",x[1]), phlist ), :entropy=>entropy_list, :log_redundancy=>redund_list, :K_complexity=>k_complex_list,
      :length_nonzeros=>len_nonzeros_list )
  open( csvfile, "w" ) do f
    hostname = readchomp(`hostname`)
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ", funcs )
    println(f,"# nreps: ", nreps )
    println(f,"# max_tries: ", max_tries )
    println(f,"# max_steps: ", max_steps )
    CSV.write( f, df, append=true, writeheader=true )
  end 
  df
end

# Compute multiset average of the entropy evolvability of a phenotype.
# Based on evolving goals to map to the phenotype
function multiset_phenotype_evolvability( p::Parameters, funcs::Vector{Func}, ph::Goal, nreps::Int64, max_tries::Int64, max_steps::Int64 )
  if max_tries <= nreps
    error("Please set max_tries to be larger than nreps in function multiset_phenotype_evolvability")
  end
  nphenos = 2^(2^p.numinputs)
  pcounts = zeros(Int64,nphenos)
  #circ_steps_list = pheno_evolve( p, funcs, ph, nreps, max_tries, max_steps )
  n_procs = nprocs()==1 ? 1 : nprocs()-1
  nreps_iter = Int(ceil( nreps/n_procs ) )
  println("nreps_iter: ",nreps_iter)
  circ_steps_list_list = pmap(_->pheno_evolve( p, funcs, ph, nreps_iter, max_tries, max_steps ), 1:n_procs )
  for circ_steps_list in circ_steps_list_list
    for (c,steps) in circ_steps_list
      phlist = mutate_all( c, funcs, output_outputs=true )
      for pph in phlist
        pcounts[pph[1]+1] += 1
      end
    end
  end
  denom = sum(pcounts)
  # println("denom: ",denom)
  #for i = 1:nphenos
  #  print(pcounts[i]," ")
  # end
  pcounts[ph[1]+1] = 0  # Set count for self-edges to be 0
  entropy = CGP.entropy( pcounts/denom)
  nonzeros = findall( x->x!=0, pcounts )
  (pcounts, entropy, length(nonzeros) )
end

# Find a random phenotype with a specified K complexity
function find_phenotype_with_Kcomplexity( p::Parameters, funcs::Vector{Func}, Kcomplexity::Int64, max_steps::Int64; kdict::Dict=Dict() )
  if length(kdict) == 0
    kdict = kolmogorov_complexity_dict(p,funcs)
  end
  step = 0
  ph = rand(MyInt(0):MyInt(2^2^p.numinputs-1))
  while step < max_steps && kdict[ph] != Kcomplexity
    ph = rand(MyInt(0):MyInt(2^2^p.numinputs-1))
    step += 1
  end
  if step == max_steps
    return nothing
  else
    return ph
  end
end

# Find a random phenotype with a specified minimum T complexity
function find_phenotype_with_min_Tcomplexity( p::Parameters, funcs::Vector{Func}, Tcomplexity::Int64, max_steps::Int64 )
  step = 0
  #ph = rand(MyInt(0):MyInt(2^2^p.numinputs-1))
  c = random_chromosome( p, funcs)
  CTcmplx = complexity5(c)
  println("CTcmplx: ",CTcmplx,"  output_values(c): ",output_values(c),"  kdict[output_values(c)[1]]: ",kdict[output_values(c)[1]],"  number active: ",number_active(c))
  while step < max_steps && CTcmplx < Tcomplexity
    c = random_chromosome( p, funcs)
    CTcmplx = complexity5(c)
    println("step: ",step,"  CTcmplx: ",CTcmplx,"  output_values(c): ",output_values(c),"  kdict[output_values(c)[1]]: ",kdict[output_values(c)[1]],"  number active: ",number_active(c))
    step += 1
  end
  #println("Tcmplx: ",CTcmplx)
  ph = output_values(c)[1]
  if step == max_steps
    println("WARNING: find_phenotype_with_Kcomplexity() failed with Kcomplexity: ",Kcomplexity," and step: ",step)
    return nothing
  else
    return ph
  end
end

# Find a random genotype with a specified Hamming distance
function find_phenotype_with_hamming_dist( p::Parameters, funcs::Vector{Func}, ph::Goal, h_dist::Int64, max_find_steps::Int64 )
  println("find_genotype_with_hamming_dist() h_dist: ",h_dist,"  ph: ",ph)
  step = 0
  rph = rand(MyInt(0):Ones)
  #hdist = hamming_distance(rph,ph[1],p.numinputs)
  #println("step: ",step,"  hdist: ",hdist)
  while step < max_find_steps && hamming_distance(rph,ph[1],p.numinputs) != h_dist/2^p.numinputs
    rph = rand(MyInt(0):Ones)
    #hdist = hamming_distance(rph,ph[1],p.numinputs)
    #println("rph: ",[rph],"  step: ",step,"  hdist: ",hdist)
    step += 1
  end
  if step == max_find_steps
    println("WARNING: find_genotype_with_hamming_dist() failed with ph: ",ph,"  h_dist: ",h_dist," and step: ",step)
    return nothing
  else
    println("find_genotype_with_hamming_dist() returned ",[rph],"  step: ",step)
    return rph
  end
end

# Find a random genotype with a specified K complexity
function find_genotype_with_Kcomplexity( p::Parameters, funcs::Vector{Func}, Kcomplexity::Int64, max_find_steps::Int64; kdict::Dict=Dict() )
  #println("find_genotype_with_Kcomplexity() Kcomplexity: ",Kcomplexity)
  if length(kdict) == 0
    kdict = kolmogorov_complexity_dict(p,funcs)
  end
  step = 0
  c = random_chromosome( p, funcs)
  while step < max_find_steps && kdict[output_values(c)[1]] != Kcomplexity
    c = random_chromosome( p, funcs)
    step += 1
  end
  if step == max_find_steps
    println("WARNING: find_genotype_with_Kcomplexity() failed with Kcomplexity: ",Kcomplexity," and step: ",step)
    return nothing
  else
    #println("find_genotype_with_Kcomplexity() kdict[output_values(c)[1]]: ",kdict[output_values(c)[1]])
    return c
  end
end

# Find a random genotype with a specified max T complexity
function find_genotype_with_max_Tcomplexity( p::Parameters, funcs::Vector{Func}, Tcomplexity::Int64, max_find_steps::Int64 )
  println("Tcomplexity: ",Tcomplexity)
  step = 0
  c = random_chromosome( p, funcs)
  CTcmplx = complexity5(c)
  println("Tcomplexity: ",Tcomplexity,"  CTcmplx: ",CTcmplx)
  while step < max_find_steps && ( CTcmplx > Tcomplexity || CTcmplx < Tcomplexity-1.0 )
    c = random_chromosome( p, funcs)
    CTcmplx = complexity5(c)
    println("step: ",step,"  CTcmplx: ",CTcmplx,"  output_values(c)[1]: ",output_values(c)[1])
    step += 1
  end
  if step == max_find_steps
    return nothing
  else
    println("step: ",step,"  CTcmplx: ",CTcmplx,"  output_values(c)[1]: ",output_values(c))
    return c
  end
end

# 
function run_evolve_to_Kcomplexity_mt( p::Parameters, funcs::Vector{Func}, fromComplexity::Int64, toComplexity::Int64, nreps::AbstractRange, max_ev_tries::Int64, 
      max_find_steps::Int64, max_evolve_steps::AbstractRange; csvfile::String="" )
  nreps_list = Int64[]
  max_ev_steps_list = Int64[]
  numfailures_list = Int64[]
  ttry_list_list = Vector{Int64}[]
  mean_steps_list = Float64[]
  median_steps_list = Float64[]
  for nr in nreps
    for mev_steps in max_evolve_steps
      push!(nreps_list,nr)
      push!(max_ev_steps_list,mev_steps)
      (numfailures, ttry_list, total_steps_list) = 
          #evolve_to_Kcomplexity( p, funcs, fromComplexity, toComplexity, nr, max_ev_tries, max_find_steps, mev_steps )
          evolve_to_Kcomplexity_mt( p, funcs, fromComplexity, toComplexity, nr, max_ev_tries, max_find_steps, mev_steps, useKcomplexity=true)
      push!(numfailures_list,numfailures)
      push!(ttry_list_list,ttry_list)
      push!(mean_steps_list,mean(total_steps_list))
      push!(median_steps_list,median(total_steps_list))
    end
  end
  len = length(nreps)*length(max_evolve_steps)    
  df = DataFrame( 
        :fromC=>fill(fromComplexity,len),
        :toC=>fill(toComplexity,len),
        :max_ev_tries=>fill(max_ev_tries,len),
        :nreps=>nreps_list,
        :failures=>numfailures_list,
        :max_ev_steps=>max_ev_steps_list,
        :mean_steps=>mean_steps_list,
        :median_steps=>median_steps_list,
        :try_list=>ttry_list_list
  )
  if length(csvfile) > 0
    hostname = readchomp(`hostname`)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", funcs )
      println(f,"# max_find_steps: ",max_find_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# Do multiple epochal evolutions ( calls to neutral_evolve() ) to evolve from a circuit of Kcomplexity fromKcomplexity to a goal of Kcomplexity toKcomplexity
function evolve_to_Kcomplexity( p::Parameters, funcs::Vector{Func}, fromKcomplexity::Int64, toKcomplexity::Int64, nreps::Int64, max_ev_tries::Int64, max_find_steps::Int64, max_evolve_steps::Int64 )
  kdict = kolmogorov_complexity_dict(p,funcs)
  ttry_list = zeros(Int64,max_ev_tries)
  total_steps_list = zeros(Int64,nreps)
  numfailures = 0
  for i = 1:nreps
    c = find_genotype_with_Kcomplexity( p, funcs, fromKcomplexity, max_find_steps, kdict=kdict )
    ph = find_phenotype_with_Kcomplexity( p, funcs, toKcomplexity, max_find_steps, kdict=kdict )
    #println("i: ",i,"  Kfrom: ",kdict[output_values(c)[1]],"  Kto: ",kdict[ph])
    @assert fromKcomplexity == kdict[output_values(c)[1]]
    @assert toKcomplexity == kdict[ph]
    ttry = 1
    (nc,steps) = neutral_evolution( c, funcs, [ph], max_evolve_steps )
    total_steps = steps
    while ttry < max_ev_tries && steps == max_evolve_steps
      (nc,steps) = neutral_evolution( c, funcs, [ph], max_evolve_steps )
      total_steps += steps
      ttry += 1
    end
    if ttry == max_ev_tries && steps == max_evolve_steps
      numfailures += 1
    end
    ttry_list[ttry] += 1
    total_steps_list[i] = total_steps
  end
  (numfailures,ttry_list,total_steps_list)
end  

# Do multiple epochal evolutions ( calls to neutral_evolve() ) to evolve from a circuit of complexity fromComplexity to a goal of complexity toComplexity
# Multi-threaded version
function evolve_to_Kcomplexity_mt( p::Parameters, funcs::Vector{Func}, fromComplexity::Int64, toComplexity::Int64, nreps::Int64, max_ev_tries::Int64, max_find_steps::Int64, max_evolve_steps::Int64;
    useKcomplexity::Bool=true ) # if useKcomplexity==false, this means use Tcomplexity
  kdict = useKcomplexity ? kolmogorov_complexity_dict(p,funcs) : Dict()
  ttry_list = [ Threads.Atomic{Int64}(0) for i= 1:max_ev_tries]
  total_steps_list = [ Threads.Atomic{Int64}(0) for i= 1:nreps]
  numfailures = Threads.Atomic{Int64}(0)
  Threads.@threads for i = 1:nreps
    if useKcomplexity
      c = find_genotype_with_Kcomplexity( p, funcs, fromComplexity, max_find_steps, kdict=kdict )
      ph = find_phenotype_with_Kcomplexity( p, funcs, toComplexity, max_find_steps, kdict=kdict )
      #println("i: ",i,"  ph: ",ph,"  output_values(c): ",output_values(c),"  kdict[output_values(c)[1]]: ",kdict[output_values(c)[1]])
      #println("i: ",i,"  Kfrom: ",kdict[output_values(c)[1]],"  Kto: ",kdict[ph])
      @assert fromComplexity == kdict[output_values(c)[1]]
      @assert toComplexity == kdict[ph]
    else # use Tononi complexity
      c = find_genotype_with_max_Tcomplexity( p, funcs, fromComplexity, max_find_steps )
      ph = find_phenotype_with_min_Tcomplexity( p, funcs, toComplexity, max_find_steps )
    end
    ttry = 1
    (nc,steps) = neutral_evolution( c, funcs, [ph], max_evolve_steps )
    total_steps = steps
    while ttry < max_ev_tries && steps == max_evolve_steps
      (nc,steps) = neutral_evolution( c, funcs, [ph], max_evolve_steps )
      total_steps += steps
      ttry += 1
    end
    if ttry == max_ev_tries && steps == max_evolve_steps
      Threads.atomic_add!( numfailures, 1 )
    end
    Threads.atomic_add!( ttry_list[ttry], 1 )
    Threads.atomic_add!( total_steps_list[i], total_steps )
  end
  (numfailures[],map(tt->tt[],ttry_list),map(ts->ts[],total_steps_list))
end  

# Returns list of unique phenotypes that are in the evolvability set of ph
function evolvability_list_evolution( p::Parameters, funcs::Vector{Func}, phlist::GoalList, ncircuits::Int64, max_tries::Int64, max_steps::Int64 )
  pheno_set = Set(Goal[])
  for ph in phlist
    circuits_steps_list = pheno_evolve( p, funcs, ph, ncircuits, max_tries, max_steps )
    circuits_list = map(x->x[1], circuits_steps_list )
    for circ in circuits_list
      phenos = mutate_all( circ, funcs, output_outputs=true, output_circuits=false )
      for phm in phenos
        push!(pheno_set,phm)
      end
    end
  end
  setdiff!( pheno_set, Set([ph for ph in phlist ]) )
  collect(pheno_set)
end

#  fit_dict = Dict{MyInt,Float64}( [(target[1],1.0)] )
function evolvability_list_filtered( p::Parameters, funcs::Vector{Func}, ph::Goal, target::Goal, evo_list::Vector{Goal}, fit_dict::Dict{MyInt,Float64}=Dict{MyInt,Float64}(); 
      neutral::Bool=false )
  if length(fit_dict) == 0
    fit_dict = Dict{MyInt,Float64}( [(target[1],1.0)] )
  end   
  ph_fit = get!( fit_dict, ph[1], rand() )
  if neutral 
    new_list = filter( x->get!( fit_dict, x[1], rand()) >= ph_fit , evo_list )
  else
    new_list = filter( x->get!( fit_dict, x[1], rand()) > ph_fit , evo_list )
  end
  if target in new_list
    println("target ",target," found")
  end
  new_list
end
