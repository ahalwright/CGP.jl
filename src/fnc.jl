using JLD, HDF5, Tables, Base.Threads

# Combines calls to enumerate_circuits(), find_neutral_components(), dict_to_csv() and consolidate_df(), 
#    and writes the resulting dataframe to a file if the csvfile keyword argument is a non-empty string.G
function component_properties( p::Parameters, pheno::MyInt, funcs::Vector{Func}=default_funcs(p.numinputs); csvfile::String="" )
  ecl = enumerate_circuits( p, funcs )
  S=find_neutral_components(ecl,pheno,funcs)
  df = dict_to_csv(S,p,funcs)
  rdf = consolidate_df(df,p,funcs,csvfile=csvfile)
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      numprocs = nprocs()==1 ? 1 : nprocs()-1
      println(f,"# host: ",hostname," with ",numprocs,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ",funcs)
      println(f,"# phenotype: ",@sprintf("0x%04x",pheno))
      println(f,"# nthreads: ",nthreads())
      CSV.write( f, rdf, append=true, writeheader=true )
    end
  end
  rdf
end

# p = Parameters(3,1,4,5)  # Example
#  funcs = default_funcs(p.numinputs) #  ecl = enumerate_circuits( p, funcs); length(ecl) 
#  @time S=find_neutral_components(ec2,0x005a); print_lengths(S)
function find_neutral_components( ch_list::Vector{Chromosome}, phenotype::MyInt, funcs::Vector{Func}=default_funcs(p.numinputs); jld_file::String="" )
  p = ch_list[1].params
  ch_list = filter( x->output_values(x)[1]==phenotype, ch_list )
  println("length(ch_list): ",length(ch_list))
  if length(ch_list) == 0
    error("no genotypes that map to the given phenothype ",[phenotype])
  end
  if nprocs() == 1
    S = find_neutral_comps( ch_list, phenotype )
  else
    # Break ch_list into sublists for parallel processing
    split = div(length(ch_list),nprocs()-1) + 1
    chl = Vector{Chromosome}[]
    for i = 0:nprocs()-2
      #println("i: ",i,"  cl: ",ch_list[i*split+1:min((i+1)*split,length(ch_list))])
      push!(chl,ch_list[i*split+1:min((i+1)*split,length(ch_list))])
    end
    #println("chl: ",chl)
    Slist = map(cl->find_neutral_comps( cl, phenotype ), chl )
    S = Dict{Int64,Set{Int128}}()
    for i = 1:length(Slist)
      S = merge_dictionaries!( S, Slist[i] )
      #println("md: S: ",S)
    end
  end
  if length(jld_file) > 0
    D = dict_int_keys_to_string_keys( S )
    save(jld_file,"D",D)
  end
  S
end

function find_neutral_comps( ch_list::Vector{Chromosome}, phenotype::MyInt )
  S = Dict{Int64,Set{Int128}}()
  new_key = 1
  for g in ch_list
    ig = chromosome_to_int(g,funcs)
    mlist = filter( x->output_values(x)[1]==phenotype, mutate_all( g, funcs, output_chromosomes=true, output_outputs=false ) )
    ihlist = map(h->chromosome_to_int(h,funcs),mlist)
    push!(ihlist,ig)
    ihset = Set(ihlist)
    if length(ihset) > 0
      for ky in keys(S)
        #println("ky: ",ky,"  S[ky]: ",S[ky])
        if length( intersect( ihset, S[ky] ) ) > 0
          union!( ihset, S[ky] )
          delete!( S, ky )
        end
      end
      S[ new_key ] = ihset
      if new_key % 100 == 0
        print("ig: ",ig,"  length(ihset): ",length(ihset),"   ")
        println("length(S[",new_key,"]) = ",length(S[new_key]))
      end
      new_key += 1
    end
  end
  for ky0 in keys(S)
    for ky1 in keys(S)
      if length(intersect(S[ky0],S[ky1]))>0 && ky0 != ky1
        println("the intersection of set S[",ky0,"] and set S[",ky1,"] is nonempty")
      end
    end
  end
  S
end

function merge_set_to_dict( set::Set{Int128}, S::Dict{Int64,Set{Int128}} )
  new_key = length(keys(S))==0 ? 1 : maximum(keys(S)) + 1
  for ky in keys(S)
    if length( intersect( set, S[ky] ) ) > 0
      union!( set, S[ky] )
      delete!( S, ky )
    end
  end
  S[ new_key ] = set
  new_key += 1
  S
end

function merge_dictionaries!( S::Dict{Int64,Set{Int128}}, S2::Dict{Int64,Set{Int128}})
  for ky in keys(S2)
    S = merge_set_to_dict( S2[ky], S )
    #println("msd: S: ",S)
  end
  S
end

function print_lengths(S)
  for ky in keys(S) 
    println("ky: ",ky,"  length(S[ky])): ",length(S[ky])) 
  end
end

function dict_int_keys_to_string_keys( S::Dict{Int64,Set{Int128}} )
  D = Dict{String,Set{Int128}}()
  for k in keys(S)
    D[@sprintf("%s",k)] = S[k]
  end
  D
end

# A less efficient but simpler test function for dict_to_csv()
function dict_csv( S::Dict{Int64,Set{Int128}}, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs) )
  key_list = Int64[]
  length_list = Int64[]
  avg_robust_list = Float64[]
  std_robust_list = Float64[]
  rng_robust_list = Float64[]
  avg_evo_list = Float64[]
  rng_evo_list = Float64[]
  std_evo_list = Float64[]
  avg_cmplx_list = Float64[]
  std_cmplx_list = Float64[]
  rng_cmplx_list = Float64[]
  for ky in keys(S)
    push!(key_list,ky)
    push!(length_list,length(S[ky]))
    robust_list = Float64[]
    evo_list = Float64[]
    cmplx_list = Float64[]
    n = length(S[ky])
    for s in S[ky]
      c = int_to_chromosome( s, p, funcs )
      (robust,evo) = mutate_all( c, funcs, robustness_only=true ) 
      cmplx = complexity5(c)
      push!(robust_list,robust)
      push!(evo_list,evo)
      push!(cmplx_list,cmplx)
    end
    push!(avg_robust_list,sum(robust_list)/n)
    push!(avg_evo_list,sum(evo_list)/n)
    push!(avg_cmplx_list,sum(cmplx_list)/n)
    push!(std_robust_list,std(robust_list))
    push!(std_evo_list,std(evo_list))
    push!(std_cmplx_list,std(cmplx_list))
    push!(rng_robust_list,maximum(robust_list)-minimum(robust_list))
    push!(rng_evo_list,maximum(evo_list)-minimum(evo_list))
    push!(rng_cmplx_list,maximum(cmplx_list)-minimum(cmplx_list))
  end
  df = DataFrame()
  df.key = key_list
  df.length = length_list
  df.avg_robust = avg_robust_list
  df.std_robust = std_robust_list
  df.rng_robust = rng_robust_list
  df.avg_evo = avg_evo_list
  df.std_evo = std_evo_list
  df.rng_evo = rng_evo_list
  df.avg_cmplx = avg_cmplx_list
  df.std_cmplx = std_cmplx_list
  df.rng_cmplx = rng_cmplx_list
  df
end

function dict_to_csv( S::Dict{Int64,Set{Int128}}, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs) )
  key_list = Int64[]
  length_list = Int64[]
  avg_robust_list = Float64[]
  std_robust_list = Float64[]
  rng_robust_list = Float64[]
  avg_evo_list = Float64[]
  rng_evo_list = Float64[]
  std_evo_list = Float64[]
  avg_cmplx_list = Float64[]
  std_cmplx_list = Float64[]
  rng_cmplx_list = Float64[]
  evo_list = Float64[]
  cmplx_list = Float64[]
  for ky in keys(S)
    push!(key_list,ky)
    push!(length_list,length(S[ky]))
    (avg_robust, std_robust, rng_robust, avg_evo, std_evo, rng_evo, avg_cmplx, std_cmplx, rng_cmplx) = robust_evo_cmplx( S[ky], p, funcs )
    push!(avg_robust_list,avg_robust)
    push!(std_robust_list,std_robust)
    push!(rng_robust_list,rng_robust)
    push!(avg_evo_list,avg_evo)
    push!(std_evo_list,std_evo)
    push!(rng_evo_list,rng_evo)
    push!(avg_cmplx_list,avg_cmplx)
    push!(std_cmplx_list,std_cmplx)
    push!(rng_cmplx_list,rng_cmplx)
  end
  df = DataFrame()
  df.key = key_list
  df.length = length_list
  df.avg_robust = avg_robust_list
  df.std_robust = std_robust_list
  df.rng_robust = rng_robust_list
  df.avg_evo = avg_evo_list
  df.std_evo = std_evo_list
  df.rng_evo = rng_evo_list
  df.avg_cmplx = avg_cmplx_list
  df.std_cmplx = std_cmplx_list
  df.rng_cmplx = rng_cmplx_list
  df
end  

function robust_evo_cmplx( set::Set{Int128}, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs) )
  sum_robust = 0.0
  sum_cmplx = 0.0
  sum_evo = 0.0
  sum_sq_robust = 0.0
  sum_sq_cmplx = 0.0
  sum_sq_evo = 0.0
  max_robust = 0.0
  max_cmplx = 0.0
  max_evo = 0.0
  min_robust = 1E9
  min_cmplx = 1E9
  min_evo = 1E9
  n = length(set)
  for s in set
    c = int_to_chromosome( s, p, funcs )
    (robust,evo) = mutate_all( c, funcs, robustness_only=true ) 
    cmplx = complexity5(c)
    sum_robust += robust; sum_sq_robust += robust^2; max_robust=robust>max_robust ? robust : max_robust; min_robust=robust<min_robust ? robust : min_robust
    sum_evo += evo; sum_sq_evo += evo^2; max_evo= evo>max_evo ? evo : max_evo; min_evo=evo<min_evo ? evo : min_evo
    sum_cmplx += cmplx; sum_sq_cmplx += cmplx^2; max_cmplx=cmplx>max_cmplx ? cmplx : max_cmplx; min_cmplx=cmplx<min_cmplx ? cmplx : min_cmplx
  end
  sd_robust = sqrt( 1.0/(n-1.0)*(sum_sq_robust - sum_robust^2/n) )
  sd_evo = sqrt( 1.0/(n-1.0)*(sum_sq_evo - sum_evo^2/n) )
  sd_cmplx = sqrt( 1.0/(n-1.0)*(sum_sq_cmplx - sum_cmplx^2/n) )
  rng_robust = max_robust-min_robust
  rng_evo = max_evo-min_evo
  rng_cmplx = max_cmplx-min_cmplx
  ( sum_robust/n, sd_robust, rng_robust, sum_evo/n, sd_evo, rng_evo, sum_cmplx/n, sd_cmplx, rng_cmplx )
end

# Consolidates df by averaging dataframe rows that correspond to rows of the same length
function consolidate_df( df::DataFrame, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs); csvfile::String="" )
  #ssum = zeros(Float64,11)
  df_matrix = Tables.matrix(df[:,3:end])
  dict = Dict{Int64,Vector{Float64}}()
  default_value = fill(-1.0,11)
  for i = 1:size(df)[1]
    len = df[i,2]
    row = vcat(Float64(len),1.0,df_matrix[i,:])  # Second element accumulates the count for this len
    prevval = get(dict,len,default_value)
    if prevval != default_value
      newval = prevval + row
      newval[1] = len  # Don't sum length column
      dict[len] = newval
      #println("len: ",len," ssum: ",newval)
    else
      dict[len] = row
      #println("len: ",len," row: ",row)
    end
  end
  ordered_keys = sort([ky for ky in keys(dict)])
  println("ordered_keys: ",ordered_keys)
  rdf = DataFrame( vcat([:len=>Float64[],:count=>Float64[]],[ Symbol(nm)=>Float64[] for nm in names(df)[3:end] ] ))
  println("names(rdf): ",names(rdf))
  for ky in ordered_keys
    #println("ky: ",ky,"  dict[ky]: ",dict[ky])
    ssum = dict[ky]
    ssum[3:end] = ssum[3:end]/ssum[2]  # Convert from sum to average
    push!(rdf,ssum)
  end
  rdf
end 

function filter_parallel( ch_list::Vector{Int64}, nprcs::Integer )
  #if nprocs() == 1
  #  filter( x->output_values(x)[1]==phenotype, ch_list )
  #else
    # Break ch_list into sublists for parallel processing
    split = div(length(ch_list),nprcs-1)
    println("split: ",split)
    chl = Vector{Int64}[]
    for i = 0:nprcs-1
      println("i: ",i,"  cl: ",ch_list[i*split+1:min((i+1)*split,length(ch_list))])
      push!(chl,ch_list[i*split+1:min((i+1)*split,length(ch_list))])
    end
    println("chl: ",chl)
    #Slist = map(cl->find_neutral_comps( cl, phenotype ), chl )
    chl
  #end
end
