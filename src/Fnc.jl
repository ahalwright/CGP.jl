using JLD, HDF5, Tables, Base.Threads

#  As of 2/20/22, julia -p 4 -L CGP.jl -L Fnc.jl -L num_active_lc.jl -L to_sublists.jl 
#  Example:  if use_lincircuit==false, let p = Parameters(3,1,4,3)  else let p = Parameters(3,1,3,2) 
#    For use_lincircuits==true, larger values of numinteriors=numinstructions and numlevelsback=numregisters exceed memory even on surt2
#  phenotype = 0x0015   # count=6848 from data/12_26_21/pheno_counts_12_26_21_F.csv (Cartesian)
#  funcs = default_funcs(p);  
#  for Chromosomes:  ecl = enumerate_circuits_ch( p, funcs); length(ecl) 
#  for LinCircuits:  ecl = enumerate_circuits_lc( p, funcs); length(ecl) 
#  phl = [0x0015,0x005a]
#  df = component_properties( pp, phl, use_lincircuit=false );
function component_properties( p::Parameters, pheno_list::Vector{MyInt}, 
      nwalks_per_set::Int64=1, walk_length::Int64=5, nwalks_per_circuit::Int64=10,
      funcs::Vector{Func}=default_funcs(p.numinputs); use_lincircuit::Bool=false, csvfile::String="", jld_file::String="" )
  function chp_list_to_rdf( chp_list, p::Parameters )
    S = find_neutral_comps( chp_list, p, funcs )
    df = dict_to_csv(S,p,chp_list[1][2],funcs,use_lincircuit=use_lincircuit,nwalks_per_set=nwalks_per_set,walk_length=walk_length,nwalks_per_circuit=nwalks_per_circuit) 
    rdf = consolidate_df(df,p,funcs)
  end
  println("walk params: ",(nwalks_per_set,walk_length,nwalks_per_circuit))
  sort!(pheno_list)
  rdf_list = DataFrame[]
  ecl = use_lincircuit ? enumerate_circuits_lc( p, funcs ) : enumerate_circuits_ch( p, funcs )
  println("length(ecl): ",length(ecl))
  chp_lists = pairs_to_sublists( ecl, pheno_list, funcs ) 
  chp_nonempty_lists = use_lincircuit ? Vector{Tuple{LinCircuit,MyInt}}[] : Vector{Tuple{Chromosome,MyInt}}[]
  for chp_list in chp_lists
    println("length(chp_list): ",length(chp_list))
    if length(chp_list) > 0
      push!(chp_nonempty_lists,chp_list)
    end
  end  
  rdf_list = pmap( chp_list->chp_list_to_rdf( chp_list, p ), chp_nonempty_lists )
  #rdf_list = map( chp_list->chp_list_to_rdf( chp_list, p ), chp_nonempty_lists )
  cdf_list = DataFrame[]
  for rdf in rdf_list
    if size(rdf)[1] >= 4
      ccdf = scorrelations(rdf)
      push!(cdf_list,ccdf)
    end
  end  
  cdf = vcat(cdf_list...)
  df = vcat(rdf_list...)
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      numprocs = nprocs()==1 ? 1 : nprocs()-1
      println(f,"# host: ",hostname," with ",numprocs,"  processes: " )
      println(f,"# funcs: ",funcs)
      println(f,"# use_lincircuit: ",use_lincircuit)
      CSV.write( f, df, append=true, writeheader=true )
      #CSV.write( f, cdf, append=true, writeheader=true )
    end
  end
  (df,cdf)
end

function find_neutral_comps( chp_list::Union{Vector{Tuple{Chromosome,MyInt}},Vector{Tuple{LinCircuit,MyInt}}}, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs) )
  #D println("find_neutral_comps: chp_list: ",chp_list)
  use_lincircuit = (typeof(chp_list)==Vector{Tuple{LinCircuit,MyInt}}) 
  S = Dict{Int64,Set{Int128}}()
  new_key = 1
  for chp in chp_list
    mlist = filter( x->output_values(x,funcs)[1]==chp[2], mutate_all( chp[1], funcs, output_circuits=true, output_outputs=false ) )
    ihlist = use_lincircuit ?  map(h->circuit_to_circuit_int(h,funcs),mlist) : map(h->chromosome_to_int(h,funcs),mlist) 
    ich = use_lincircuit ? circuit_to_circuit_int(chp[1],funcs) : chromosome_to_int(chp[1],funcs)
    push!(ihlist,ich)
    ihset = Set(ihlist)
    #D println("ich: ",ich,"  ihset: ",ihset)
    ihphenos = use_lincircuit ? map(ic->output_values(circuit_int_to_circuit(ic, p, funcs))[1],ihlist) : map(ic->output_values(int_to_chromosome(ic,p,funcs))[1],ihlist) 
    #D println("ihphenos: ",ihphenos)   # Correctness check
    if length(ihset) > 0
      for ky in keys(S)
        ##D println("ky: ",ky,"  S[ky]: ",S[ky])
        if length( intersect( ihset, S[ky] ) ) > 0
          union!( ihset, S[ky] )
          ##D println("ihset after union!: ",ihset)
          delete!( S, ky )
        end
      end
      S[ new_key ] = ihset
      if new_key % 100 == 0
        #D print("ich: ",ich,"  length(ihset): ",length(ihset),"   ")
        #D println("length(S[",new_key,"]) = ",length(S[new_key]))
      end
      new_key += 1
    end
  end
  for ky0 in keys(S)
    for ky1 in keys(S)
      if length(intersect(S[ky0],S[ky1]))>0 && ky0 != ky1
        #D println("the intersection of set S[",ky0,"] and set S[",ky1,"] is nonempty")
      end
    end
  end
  S
end

function dataframe_count_phenos( rdf::DataFrame )
  @assert "len" in names(rdf)
  @assert "count" in names(rdf)
  ssum = 0
  for i = 1:size(rdf)[1]
    ssum += rdf.len[i]*rdf.count[i]
  end
  ssum
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

function dict_to_csv( S::Dict{Int64,Set{Int128}}, p::Parameters, pheno::MyInt, funcs::Vector{Func}=default_funcs(p.numinputs); 
    use_lincircuit::Bool=false, nwalks_per_set::Int64=20, walk_length::Int64=50, nwalks_per_circuit::Int64=3 )
  #println("dict_to_csv nwalks_per_set: ",nwalks_per_set,"  walk_length: ",walk_length,"  nwalks_per_circuit: ",nwalks_per_circuit)
  #println("dict_to_csv: use_lincircuit: ",use_lincircuit)
  #println("S: ",S)
  key_list = Int64[]
  length_list = Int64[]
  avg_robust_list = Float64[]
  std_robust_list = Float64[]
  rng_robust_list = Float64[]
  avg_evo_list = Float64[]
  rng_evo_list = Float64[]
  std_evo_list = Float64[]
  if !use_lincircuit 
    cmplx_list = Float64[]
    avg_cmplx_list = Float64[]
    std_cmplx_list = Float64[]
    rng_cmplx_list = Float64[]
  end
  evo_list = Float64[]
  avg_walk_list = Float64[]
  sum_ma_walk_list = Float64[]
  avg_numactive_list = Float64[]
  walk_count = 1
  for ky in keys(S)
    push!(key_list,ky)
    push!(length_list,length(S[ky]))
    if !use_lincircuit 
      (avg_robust, std_robust, rng_robust, avg_evo, std_evo, rng_evo, avg_cmplx, std_cmplx, rng_cmplx, avg_walk, sum_ma_walk, avg_numactive ) = 
        robust_evo_cmplx( S[ky], p, funcs, use_lincircuit=use_lincircuit, nwalks_per_set=nwalks_per_set, walk_length=walk_length, nwalks_per_circuit=nwalks_per_circuit )
    else
      (avg_robust, std_robust, rng_robust, avg_evo, std_evo, rng_evo, avg_walk, sum_ma_walk, avg_numactive ) = 
        robust_evo_cmplx( S[ky], p, funcs, use_lincircuit=use_lincircuit, nwalks_per_set=nwalks_per_set, walk_length=walk_length, nwalks_per_circuit=nwalks_per_circuit )
    end
    #println("dict_to_csv: walk_count: ",walk_count,"   avg_walk: ",avg_walk)
    push!(avg_robust_list,avg_robust)
    push!(std_robust_list,std_robust)
    push!(rng_robust_list,rng_robust)
    push!(avg_evo_list,avg_evo)
    push!(std_evo_list,std_evo)
    push!(rng_evo_list,rng_evo)
    if !use_lincircuit 
      push!(avg_cmplx_list,avg_cmplx)
      push!(std_cmplx_list,std_cmplx)
      push!(rng_cmplx_list,rng_cmplx)
    end
    push!(avg_walk_list, avg_walk)
    push!(sum_ma_walk_list,sum_ma_walk)
    push!(avg_numactive_list,avg_numactive)
    walk_count += 1
  end
  lendf = length(key_list)
  df = DataFrame()
  df.key = key_list
  df.length = length_list
  df.pheno = fill(pheno,lendf)
  df.numinputs = fill(p.numinputs,lendf)
  if use_lincircuit
    df.ninstr = fill(p.numinteriors,lendf)
    df.nregs = fill(p.numlevelsback,lendf)
  else
    df.ngates = fill(p.numinteriors,lendf)
    df.levsback = fill(p.numlevelsback,lendf)
  end
  df.nwalks_set = fill(nwalks_per_set,lendf)
  df.walk_length = fill(walk_length,lendf)
  df.nwalks_circ = fill(nwalks_per_circuit,lendf)
  #df.pheno = fill(pheno,lendf)
  df.avg_robust = avg_robust_list
  #df.std_robust = std_robust_list
  #df.rng_robust = rng_robust_list
  df.avg_evo = avg_evo_list
  #df.std_evo = std_evo_list
  #df.rng_evo = rng_evo_list
  if !use_lincircuit 
    df.avg_cmplx = avg_cmplx_list
    #df.std_cmplx = std_cmplx_list
    #df.rng_cmplx = rng_cmplx_list
  end
  df.avg_walk = avg_walk_list
  df.sum_ma_walk = sum_ma_walk_list
  df.avg_nactive = avg_numactive_list
  df
end   # dict_to_csv

function robust_evo_cmplx( set::Set{Int128}, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs); 
    use_lincircuit::Bool=false, nwalks_per_set::Int64=30, walk_length::Int64=50, nwalks_per_circuit::Int64=3 )
  #println("robust_evo_cmplx use_lincircuit: ",use_lincircuit,"   nwalks_per_set: ",nwalks_per_set,"  walk_length: ",walk_length)
  sum_robust = 0.0
  sum_cmplx = 0.0
  sum_evo = 0.0
  sum_sq_robust = 0.0
  if !use_lincircuit sum_sq_cmplx = 0.0 end
  sum_sq_evo = 0.0
  max_robust = 0.0
  if !use_lincircuit max_cmplx = 0.0 end
  max_evo = 0.0
  min_robust = 1E9
  if !use_lincircuit min_cmplx = 1E9 end
  min_evo = 1E9
  sum_walk = 0
  sum_ma_walk = 0 
  n = length(set)
  walk_count = 1
  sum_numactive = 0
  for s in set
    if use_lincircuit
      c = circuit_int_to_circuit( s, p, funcs )
    else
      c = int_to_chromosome( s, p, funcs )
    end     
    (robust,evo) = mutate_all( c, funcs, robustness_only=true ) 
    cmplx = use_lincircuit ? 0 : complexity5(c)
    sum_robust += robust; sum_sq_robust += robust^2; max_robust=robust>max_robust ? robust : max_robust; min_robust=robust<min_robust ? robust : min_robust
    sum_evo += evo; sum_sq_evo += evo^2; max_evo= evo>max_evo ? evo : max_evo; min_evo=evo<min_evo ? evo : min_evo
    if !use_lincircuit 
      sum_cmplx += cmplx; sum_sq_cmplx += cmplx^2; max_cmplx=cmplx>max_cmplx ? cmplx : max_cmplx; min_cmplx=cmplx<min_cmplx ? cmplx : min_cmplx
    end
    sum_walk_plus = walk_count <= nwalks_per_set ? rand_evo_walks(deepcopy(c),walk_length,nwalks_per_circuit,funcs)[end]/nwalks_per_circuit : 0
    sum_walk += sum_walk_plus
    sum_ma_walk_plus = walk_count <= nwalks_per_set ? rand_evo_walks_mutate_all(deepcopy(c),walk_length,nwalks_per_circuit,funcs)[end]/nwalks_per_circuit : 0
    sum_ma_walk += sum_ma_walk_plus
    walk_count += 1
    numactive = use_lincircuit ? num_active_lc( c, funcs ) : number_active_gates( c )
    sum_numactive += numactive
  end
  sd_robust = sqrtn( 1.0/(n-1.0)*(sum_sq_robust - sum_robust^2/n) )
  sd_evo = sqrtn( 1.0/(n-1.0)*(sum_sq_evo - sum_evo^2/n) )
  if !use_lincircuit sd_cmplx = sqrtn( 1.0/(n-1.0)*(sum_sq_cmplx - sum_cmplx^2/n) ) end
  rng_robust = max_robust-min_robust
  rng_evo = max_evo-min_evo
  if !use_lincircuit rng_cmplx = max_cmplx-min_cmplx end
  #( sum_robust/n, sd_robust, rng_robust, sum_evo/n, sd_evo, rng_evo, sum_cmplx/n, sd_cmplx, rng_cmplx, sum_walk, sum_mutall_walk )
  #println("robust_evo_complx(): sum_walk: ",sum_walk,"  avg_walk: ",sum_walk/min(n,nwalks_per_set))
  if !use_lincircuit 
    ( sum_robust/n, sd_robust, rng_robust, sum_evo/n, sd_evo, rng_evo, sum_cmplx/n, sd_cmplx, rng_cmplx, 
        sum_walk/min(n,nwalks_per_set), sum_ma_walk/min(n,nwalks_per_set), sum_numactive/n )
  else
    ( sum_robust/n, sd_robust, rng_robust, sum_evo/n, sd_evo, rng_evo, sum_walk/min(n,nwalks_per_set), sum_ma_walk/min(n,nwalks_per_set), sum_numactive/n )
  end
end

# square root function that returns 0 for negative numbers.  Prevents errors due to roundoff.
sqrtn( x::Float64 ) = x >= 0.0 ? sqrt(x) : 0.0

# Consolidates df by averaging dataframe rows that correspond to rows of the same value of len
# Also, converts the pheno column which is Float64 to string by using the prhex() function
function consolidate_df( df::DataFrame, p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs) )
  ssum = zeros(Float64,size(df)[2]-1)
  df_matrix = Tables.matrix(df[:,3:end])
  dict = Dict{Int64,Vector{Float64}}()
  default_value = fill(-1.0,size(df)[2]-1)
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
  #println("names(rdf): ",names(rdf))
  for ky in ordered_keys
    #println("ky: ",ky,"  dict[ky]: ",dict[ky])
    ssum = dict[ky]
    ssum[3:end] = ssum[3:end]/ssum[2]  # Convert from sum to average
    #println("consolidate: length(ssum): ",length(ssum),"  size(rdf): ",size(rdf))  # uncomment if there is a row length conflict
    push!(rdf,ssum)
  end
  rdf.pheno = map(x->prhex(MyInt(Int(x))),rdf.pheno)
  rdf
end 

# Never called and not fully debugged
function filter_parallel( ch_list::Vector{Int64}, nprcs::Integer )
  #if nprocs() == 1
  #  filter( x->output_values(x,funcs)[1]==phenotype, ch_list )
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

# replaced by pheno_counts_ch() on about 12/28/21
function pheno_counts( p::Parameters, funcs::Vector{Func}; csvfile::String="" )
  pheno_counts_ch( p, funcs, csvfile=csvfile )
end

# Returns a DataFrame with 2 columns: goal and counts.  
#   counts is the exact number of genotypes which map to the goal phenotype.
# If output_vect==true, returns a pair where the first element of the pair is the DataFrame described above, 
#   and the second element of the pair is a vector P so that P(i) is the phenotype mapped to by the
#    circuit (genotype) corresponding to chromosome_int i.
function pheno_counts_ch( p::Parameters, funcs::Vector{Func}; csvfile::String="", output_vect::Bool=false )
  counts = zeros(Int64,2^2^p.numinputs)
  eci = collect(0:Int64(ceil(count_circuits_ch( p, nfuncs=length(funcs))) ) )
  if output_vect 
    P = fill(MyInt(0),length(eci))
  end
  for ich in eci
    ch = int_to_chromosome( ich, p, funcs )
    indx = output_values(ch,funcs)[1]
    counts[indx+1] = counts[indx+1] + 1
    if output_vect 
      P[ich+1] = indx
    end
  end
  df = DataFrame()
  df.goal = map(g->@sprintf("0x%04x",g),map(MyInt,collect(0:2^2^p.numinputs-1)))
  df.counts = counts
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      numprocs = nprocs()==1 ? 1 : nprocs()-1
      println(f,"# host: ",hostname," with ",numprocs,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ",funcs)
      println(f,"# nthreads: ",nthreads())
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  return output_vect ?  (df,P) :  df
end

# Returns a DataFrame with 2 columns: goal and counts.  
#   counts is the exact number of genotypes which map to the goal phenotype.
# If output_vect==true, returns a pair where the first element of the pair is the DataFrame described above, 
#   and the second element of the pair is a vector P so that P(i) is the phenotype mapped to by the
#    circuit (genotype) corresponding to circuit_int i.
function pheno_counts_lc( p::Parameters, funcs::Vector{Func}; csvfile::String="", output_vect::Bool=false )
  println("pheno_counts_lc  p: ",p)
  counts = zeros(Int64,2^2^p.numinputs)
  lci = collect(0:(count_circuits_lc( p, nfuncs=length(funcs))-1) )
  if output_vect 
    P = fill(MyInt(0),length(lci))
  end
  for ilc in lci
    lc = circuit_int_to_circuit( Int128(ilc), p, funcs )
    indx = output_values(lc,funcs)[1]
    #println("ilc: ",ilc,"  indx: ",@sprintf("0x%04x",indx))
    counts[indx+1] = counts[indx+1] + 1
    if output_vect 
      P[ilc+1] = indx
    end
  end
  df = DataFrame()
  df.goal = map(g->@sprintf("0x%04x",g),map(MyInt,collect(0:2^2^p.numinputs-1)))
  df.counts = counts
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      numprocs = nprocs()==1 ? 1 : nprocs()-1
      println(f,"# host: ",hostname," with ",numprocs,"  processes: " )
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ",funcs)
      println(f,"# nthreads: ",nthreads())
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  return output_vect ?  (df,P) :  df
end

# If normalize, divide each row of the matrix by the redundancy of the corresponding phenotype
function pheno_network_matrix_df( p::Parameters, funcs::Vector{Func}; normalize::Bool=false )
  nphenos = 2^(2^p.numinputs)
  phnet_matrix = zeros( Int64, nphenos, nphenos )
  (pdf, circ_to_phenotype) = pheno_counts_ch( p, funcs, output_vect=true )
  for circint = 0:length(circ_to_phenotype)-1
    circ = int_to_chromosome( circint, p, funcs )
    from_ph = output_values( circ )[1]
    mut_phenos = map(x->x[1], mutate_all( circ, funcs ) )
    for to_ph in mut_phenos
      phnet_matrix[ from_ph+MyInt(1), to_ph+MyInt(1) ] += 1
    end
  end
  if normalize 
    phnet_matrix_float = map(x->Float64(x),phnet_matrix)
    for i = 1:size(phnet_matrix_float)[1]
      phnet_matrix_float[i,:] /= pdf.counts[i]
      phnet_matrix_float[i,i] = 0.0
    end
    phdf = matrix_to_dataframe( phnet_matrix_float, goallist, hex=true, redund_column=pdf.counts )
  else
    phdf = matrix_to_dataframe( phnet_matrix, goallist, hex=true, redund_column=pdf.counts )
  end
  open( csvfile, "w" ) do f
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ", Main.CGP.default_funcs(p))
    CSV.write( f, phdf, append=true, writeheader=true )
  end
  phdf
end

function pheno_network_matrix( p::Parameters, funcs::Vector{Func} )
  nphenos = 2^(2^p.numinputs)
  phnet_matrix = zeros( Int64, nphenos, nphenos )
  (pdf, circ_to_phenotype) = pheno_counts_ch( p, funcs, output_vect=true )
  for circint = 0:length(circ_to_phenotype)-1
    circ = int_to_chromosome( circint, p, funcs )
    from_ph = output_values( circ )[1]
    mut_phenos = map(x->x[1], mutate_all( circ, funcs ) )
    for to_ph in mut_phenos
      phnet_matrix[ from_ph+MyInt(1), to_ph+MyInt(1) ] += 1
    end
  end
  phnet_matrix
end  

# Never called
function random_walk_repeats( ch::Chromosome, steps::Int64, maxsteps::Int64, funcs::Vector{Func})
  @assert ch.params.numoutputs == 1
  counts = Dict{Int128,Int64}()
  goal = output_values(ch,funcs)[1]
  println("goal: ",@sprintf("0x%04x",goal))
  mlist = mutate_all( ch, funcs )
  if length(mlist) == 0
    println("chromosome ch in random_walk_repeats() has robustness 0")
    return counts
  end
  total_steps = 0
  for i = 1:steps
    j = 1
    while j <= maxsteps   # terminated by a break statement
      sav_ch = deepcopy(ch)
      (ch,active) = mutate_chromosome!( ch, funcs )
      output = output_values(ch,funcs)[1]
      if output == goal
        #println("successful step for i= ",i,"  output: ",output)
        break
      end
      ch = sav_ch
      j += 1
    end
    if j < maxsteps
      ci = chromosome_to_int( ch, funcs )
      res = get( counts, ci, -1 )
      if res == -1
        counts[ci] = 1
      else
        counts[ci] = res + 1
      end
    else
      error("Failed to find neutral mutation")
    end
    total_steps += j
  end
  println("total_steps: ",total_steps)  
  total_counts = Dict{Int64,Int64}()
  for ky in keys(counts)
    cnt = counts[ky]
    res = get( total_counts, cnt, -1 )
    if res == -1
      total_counts[cnt] = 1
    else
      total_counts[cnt] = res + 1
    end
  end
  total_counts
end

# Do random neutral walks starting with chromosome/circuit c and return the average number of new unique genotypes per step
# Examples:
# rand_evo_walks( random_chromosome(p,funcs), 10, 2 )
# rand_evo_walks( rand_lcircuit(p,funcs), 10, 2 )
function rand_evo_walks( c::Union{Chromosome,LinCircuit}, walk_length::Int64, nwalks_per_circuit::Int64, funcs::Vector{Func} )
  #D print("rand_evo_walks ic: ",circuit_to_circuit_int(c,funcs),"   ")
  count_mutations = 0
  sum_unique_genotypes = zeros(Int64,walk_length)
  use_lincircuit = typeof(c)==LinCircuit
  phenotype = output_values(c,funcs)
  #D println("phenotype: ",prhex(phenotype[1]))
  returned_genotypes = Set{Int128}()
  new_length = 0
  new_c = deepcopy(c)
  for w = 1:nwalks_per_circuit
    #D println("walk: ",w)
    for step = 1:walk_length
      if output_values(new_c,funcs) == phenotype
        ic = use_lincircuit ? circuit_to_circuit_int(new_c,funcs) : chromosome_to_int(c,funcs) 
        union!( returned_genotypes, Set([ic]))
        #D println("step: ",step,"  returned_genotypes: ",returned_genotypes)
        prev_length = new_length
        new_length = length(returned_genotypes)
        sum_unique_genotypes[step] += (new_length - prev_length)
        #D println("step: ",step,"  new_len-prev_len: ",new_length-prev_length, "  sum_unique: ",sum_unique_genotypes[step] )
        #D println("  genotypes: ",returned_genotypes)
        c = new_c
      end
      new_c = deepcopy(c)
      use_lincircuit ? mutate_circuit!(new_c,funcs) : mutate_chromosome!(new_c,funcs)
      count_mutations += 1
    end
  end
  #println("return: ",length(returned_genotypes))
  length(returned_genotypes)
end

# Do random neutral walks starting with chromosome/circuit c and return the average number of new unique genotypes per step
# Examples:
# rand_evo_walks_mutate_all( random_chromosome(p,funcs), 10, 2 )
# rand_evo_walks_mutate_all( rand_lcircuit(p,funcs), 10, 2 )
function rand_evo_walks_mutate_all( c::Union{Chromosome,LinCircuit}, walk_length::Int64, nwalks_per_circuit::Int64, funcs::Vector{Func}) 
  count_mutations = 0
  p = c.params
  sum_unique_genotypes = zeros(Int64,walk_length)
  use_lincircuit = typeof(c)==LinCircuit
  phenotype = output_values(c,funcs)
  #D println("phenotype: ",prhex(phenotype[1]))
  for w = 1:nwalks_per_circuit
    #D println("walk: ",w)
    returned_genotypes = Set{Int128}()
    for step = 1:walk_length
      phenos,circuits = mutate_all(c,funcs,output_outputs=true,output_circuits=true)
      count_mutations += length(circuits)
      circuits = robust_filter(phenotype[1],circuits,phenos)
      #D println("step: ",step,"  circuits: ",circuits)
      #circuit_ints = use_lincircuit ? map(x->circuit_to_circuit_int(LinCircuit(x,p),funcs),circuits) : map(x->chromosome_to_int(x,funcs),circuits) 
      circuit_ints = use_lincircuit ? map(x->circuit_to_circuit_int(x,funcs),circuits) : map(x->chromosome_to_int(x,funcs),circuits) 
      union!(returned_genotypes,Set(circuit_ints))
      #D println("step: ",step,"  unique: ",length(returned_genotypes) )
      sum_unique_genotypes[step] += length(returned_genotypes)
      #D print("step: ",step,"  sum_unique: ",sum_unique_genotypes[step] )
      use_lincircuit ? mutate_circuit!(c,funcs) : mutate_chromosome!(c,funcs)
      #D println()
    end
  end
  #D println("count_mutations: ",count_mutations)
  sum_unique_genotypes./nwalks_per_circuit
end

function robust_filter( pheno::MyInt, circuits::Union{Vector{Chromosome},Vector{LinCircuit}}, phenos::Vector{Vector{MyInt}} )
  @assert length(circuits) == length(phenos)
  new_circuits = Union{Chromosome,LinCircuit}[]
  for i = 1:length(circuits)
    if phenos[i][1] == pheno
      push!(new_circuits,circuits[i])
    end
  end
  new_circuits
end

# Returns a MyInt written in hex notataion as a string
# Works for MyInt = UInt8, UInt16 but not UInt32 (TODO: fix)
prhex( x::MyInt ) = @sprintf("0x%04x",x)

function scorrelations( rdf::DataFrame )
  df = DataFrame()
  df.pheno = [rdf.pheno[1]]
  df.count = [rdf.len'*rdf.count]
  df.ncomps = [size(rdf)[1]]
  #df.numinputs = [rdf.numinputs[1]]
  #df.ngates = [rdf.ngates[1]]
  #df.levsback = [rdf.levsback[1]]
  df.walklen = [rdf.walk_length[1]]
  df.nwlkset = [rdf.nwalks_set[1]]
  df.nwlkcirc = [rdf.nwalks_circ[1]]
  df.rbst = [spearman_cor(rdf,:len,:avg_robust)[1]]
  df.rbstp = [spearman_cor(rdf,:len,:avg_robust)[2]]
  df.evo = [spearman_cor(rdf,:len,:avg_evo)[1]]
  df.evop = [spearman_cor(rdf,:len,:avg_evo)[2]]
  df.cplx = [spearman_cor(rdf,:len,:avg_cmplx)[1]]
  df.cplxp = [spearman_cor(rdf,:len,:avg_cmplx)[2]]
  df.walk = [spearman_cor(rdf,:len,:avg_walk)[1]]
  df.walkp = [spearman_cor(rdf,:len,:avg_walk)[2]]
  df.mawlk = [spearman_cor(rdf,:len,:sum_ma_walk)[1]]
  df.mawlkp = [spearman_cor(rdf,:len,:sum_ma_walk)[2]]
  df.nactive = [spearman_cor(rdf,:len,:avg_nactive)[1]]
  df.nactivep = [spearman_cor(rdf,:len,:avg_nactive)[2]]
  df
end

function scorrelations_df( rdf::DataFrame )
  current_pheno = rdf.pheno[1]
  current_pheno_index = 1
  println("  current_pheno: ",current_pheno)
  #=
  cdf = DataFrame()
  symnames = map(Symbol,names(rdf))
  for j = 1:size(rdf)[2]
    insertcols!(cdf,symnames[j]=>Vector{typeof(rdf[1,j])}[] )
  end
  =#
  ssdf = DataFrame()  # Establish scope
  sdf = DataFrame()   # Establish scope
  k = 1
  for i = 2:size(rdf)[1]
    println("i: ",i,"  rdf.pheno[i]: ",rdf.pheno[i])
    if current_pheno != rdf.pheno[i] || i == size(df)[1]
      ndf = rdf[current_pheno_index:i,:]
      if size(ndf)[1] > 4
        sdf = scorrelations( ndf )
        if k == 1
          ssdf = sdf
        else
          ssdf = vcat(ssdf,sdf)
        end
        k += 1
      end
      println("i: ",i,"  size(sdf): ",size(sdf),"  size(ssdf): ",size(ssdf))
      current_pheno = rdf.pheno[i]
      current_pheno_index = i
    end
  end
  ssdf
end

# included from to_sublists.jl on 4/10/22
using Test

function pairs_to_sublists( ecl::Union{Vector{Chromosome},Vector{LinCircuit}}, pheno_list::Vector{MyInt}, funcs::Vector{Func} )
  sort!(pheno_list)
  chp_list = [ (ec,output_values(ec,funcs)[1]) for ec in ecl ] 
  sort!( chp_list, by=x->x[2] )  
  if typeof(ecl) == Vector{Chromosome}
    ch_pheno_lists = Vector{Tuple{Chromosome,MyInt}}[]
    empty_chp_sublist = Tuple{Chromosome,MyInt}[]
  elseif typeof(ecl) == Vector{LinCircuit}
    ch_pheno_lists =Vector{Tuple{LinCircuit,MyInt}}[]
    empty_chp_sublist = Tuple{Union{Chromosome,LinCircuit},MyInt}[]
  end
  chp_sublist = deepcopy(empty_chp_sublist)
  if length(pheno_list)==0 return ch_pheno_lists end
  #D println("pheno_list: ",pheno_list)
  i = 1
  j = 1
  chp = chp_list[j]
  while i <= length(pheno_list) && j <= length(chp_list)
    chp = chp_list[j]
    #D println("j: ",j,"  i: ",i,"  chp[2]: ",prhex(chp[2]),"  pheno_list[i]: ",prhex(pheno_list[i]),"  length(chp_sublist): ",length(chp_sublist))
    if chp[2] < pheno_list[i]
       j += 1  
       #D print("j: ",j,"  increment j  ")
    elseif chp[2] == pheno_list[i]
      push!(chp_sublist,chp)
      #D print("j: ",j,"  chp: ",prhex(chp[2]),"  push chp:  length chp_sublist:  ",length(chp_sublist),"   ")
      j += 1
    elseif chp[2] > pheno_list[i]
      push!(ch_pheno_lists,chp_sublist)                
      #D print("j: ",j,"  i: ",i,"  chp: ",prhex(chp[2]),"  push chp_sublist: length chp_sublist: ",length(chp_sublist))
      #D println("  length ch_pheno_lists: ",length(ch_pheno_lists))
      i += 1
      chp_sublist = deepcopy(empty_chp_sublist)
    end
  end  # while
  #D println()
  #D println("end while: i: ",i,"  j: ",j)
  while i <= length(pheno_list)
    push!(ch_pheno_lists,chp_sublist)                
    #D print("j: ",j,"  i: ",i,"  chp: ",prhex(chp[2]),"  push chp_sublist: length chp_sublist: ",length(chp_sublist))
    #D println("  length ch_pheno_lists: ",length(ch_pheno_lists))
    i += 1
  end
  if j > length(chp_list)
    push!(ch_pheno_lists,chp_sublist)                
    #D print("j: ",j,"  i: ",i,"  chp: ",prhex(chp[2]),"  push chp_sublist: length chp_sublist: ",length(chp_sublist))
  end
  ch_pheno_lists
end

function test_to_sublists( ecl::Union{Vector{Chromosome},Vector{LinCircuit}}, pheno_list::Vector{MyInt}, funcs::Vector{Func} )
  #@testset begin "testing to_sublists() using random sets of phenotypes"
    chp_lists=pairs_to_sublists( ecl, pheno_list, funcs )
    sublists_lengths = map(x->length(x),chp_lists)
    #D println("sublists_lengths: ",sublists_lengths," length(sublists_lengths): ",length(sublists_lengths))
    i = 1
    for ph in pheno_list
      ch_list = filter( x->output_values(x,funcs)[1]==ph, ecl )
      #D println("i: ",i,"  ph: ",prhex(ph),"  length(ch_list): ",length(ch_list))
      #@test length(ch_list) == sublists_lengths[i]
      i += 1
    end
  #end # testset
end

function random_phl( p::Parameters, prob::Float64 )
  phl = MyInt[]
  for i = 0:(2^2^p.numinputs-1)
    if rand() < prob
      push!(phl,MyInt(i))
    end
  end
  phl
end

# Returns the number active not inluding inputs.
function num_active_lc( circ::LinCircuit, funcs::Vector{Func}=default_funcs(circ.params) )
  p = circ.params
  numinstructions = p.numinteriors
  res = num_active_lc( circ, numinstructions, numinstructions+1, funcs )
  reduce(+,res[2])
end

# Included from num_active_lc.jl
function num_active_lc( circ::LinCircuit, i::Int64, count::Int64, funcs::Vector{Func}=default_funcs(circ.params) )
  if count <= 0   # return if recursion infinite loop
    return
  end
  #D println("num_active_lc: i: ",i,"  count: ",count)
  p = circ.params
  instruction_active = fill(false,p.numinteriors)  # p.numinteriors == numinstructions
  v = circ.circuit_vects
  @assert p.numoutputs == 1
  n_instructions = length(circ.circuit_vects)
  Rinit = fill(MyInt(0), p.numlevelsback + p.numinputs ) # numlevelsback is the number of computational registers
  Rinit[(p.numlevelsback+1):end] = construct_context(p.numinputs)
  if v[i][2] == 1    # instruction does output result to output register
    rval = eval_lincircuit( circ, i, count-1, instruction_active, funcs )
    #D println("v[i][2] == 1: count: ",count,"  rval: ",rval)
    return ( rval, instruction_active )
  else  # instruction does not output result to output register: search for last instruction that does
    j = i-1
    while j >= 1 && v[j][2] != 1
      j -= 1
    end
    if j >=1
      instruction_active[j] = true
      rval = eval_lincircuit( circ, j, count-1, instruction_active, funcs )
    else
      rval = Rinit[1]
    end
    #D println("j: ",j,"  rval: ",prhex(rval))
    return ( rval, instruction_active )
  end
end

function eval_lincircuit( circ::LinCircuit, i::Int64, count::Int64, instruction_active::Vector{Bool},
      funcs::Vector{Func}=default_funcs(circ.params) ) 
  #D println("eval_lincircuit: i: ",i,"  count: ",count)
  instruction_active[i] = true
  p = circ.params
  v = circ.circuit_vects
  Rinit = fill(MyInt(0), p.numlevelsback + p.numinputs ) # numlevelsback is the number of computational registers
  Rinit[(p.numlevelsback+1):end] = construct_context(p.numinputs)
  done = false   # establish scope
  rc = MyInt[]
  for m = 3:4
    #println("m: ",m)
    done = false
    if v[i][m] >= 3
      #D println("i: ",i,"  m: ",m,"  v[i][m]: ",v[i][m],"  push context Rinit[v[i][m]]: ",prhex(Rinit[v[i][m]]))
      push!(rc,Rinit[v[i][m]])
      done = true
    else
      for j = (i-1):-1:1
        if v[j][2] == v[i][m]
          #D println("i: ",i,"  before recurse: rc: ",rc)
          rval = eval_lincircuit( circ, j, count-1, instruction_active, funcs )
          #D println("i: ",i,"  j: ",j,"  rc: ",rc,"  push rval: ",prhex(rval))
          push!(rc,rval)
          done = true
          break
        #else
          #println("i: ",i," v[i][m]: ",v[i][m],"  j: ",j," push Rinit[v[i][m]]: ",prhex(Rinit[v[i][m]]))
          #push!(rc,Rinit[v[i][m]])
          #done = true
        end
      end
    end
    if !done
      #D println("i: ",i," m: ",m,"  !done push Rinit[v[i][m]]: ",prhex(Rinit[v[i][m]]))
      push!(rc,Rinit[v[i][m]])
    end
    #D println("i: ",i,"  m: ",m,"  rc: ",rc)
  end
  return_value = funcs[circ.circuit_vects[i][1]].func(rc[1],rc[2])
  #D println("i: ",i,"  func: ",funcs[circ.circuit_vects[i][1]],"  return value: ",prhex(return_value),"  count: ",count)
  return return_value
end
  
prhex( x::MyInt ) = @sprintf("0x%04x",x)
  
function execl( circuit::LinCircuit, funcs::Vector{Func} )
  p = circuit.params
  R = fill(MyInt(0), p.numlevelsback + p.numinputs ) # numlevelsback is the number of computational registers
  R[(p.numlevelsback+1):end] = construct_context(p.numinputs)
  for lc in circuit.circuit_vects
    R[lc[2]] = funcs[lc[1]].func(R[lc[3]],R[lc[4]])
    #D println("lc: ",lc,"  R: ",R)
  end
  R
end

function test_num_active_lc( p::Parameters, funcs::Vector{Func}, numtries::Int64 )
  numactive_counts = zeros( Int64, p.numinteriors+1 )
  print("test_num_active: p: ",p,"  numtries: ",numtries)
  done = false
  while numtries > 0 && !done
    numinstructions = p.numinteriors
    circ = rand_lcircuit( p, funcs )
    Rr = execute_lcircuit( circ, funcs )
    nalc = num_active_lc(circ,numinstructions,numinstructions+1,funcs)
    numactive = reduce(+,nalc[2])
    numactive_counts[ numactive+1 ] += 1
    done = Rr[1]==nalc[1] ? false : true
    numtries -= 1
  end
  println("  numtries: ",numtries)
  println("numactive_counts: ",numactive_counts')
end 
