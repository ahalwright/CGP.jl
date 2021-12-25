using JLD, HDF5
# p = Parameters(3,1,4,5)  # Example
#  funcs = default_funcs(p.numinputs) #  ec2 = enumerate_circuits( p, funcs); length(ec2) 
#  @time S=find_neutral_components(ec2,0x005a); print_lengths(S)
function find_neutral_components( ch_list::Vector{Chromosome}, phenotype::MyInt; jld_file::String="" )
  p = ch_list[1].params
  funcs = default_funcs(p.numinputs)
  ch_list = filter( x->output_values(x)[1]==phenotype, ch_list )
  println("length(ch_list): ",length(ch_list))
  if length(ch_list) == 0
    error("no genotypes that map to the given phenothype ",[phenotype])
  end
  #readline()
  if nprocs() == 1
    S = find_neutral_comps( ch_list, phenotype )
  else
    split = div(length(ch_list),nprocs()-1) + 1
    # Break ch_list into sublists for parallel processing
    chl = Vector{Chromosome}[]
    for i = 0:nprocs()-2
      #println("i: ",i,"  cl: ",ch_list[i*split+1:min((i+1)*split,length(ch_list))])
      push!(chl,ch_list[i*split+1:min((i+1)*split,length(ch_list))])
    end
    #println("chl: ",chl)
    Slist = map(cl->find_neutral_comps( cl, phenotype ), chl )
    for i = 1:length(Slist)
      println("Slist[",i,"]: ",Slist[i])
    end
    S = Dict{Int64,Set{Int128}}()
    for i = 1:length(Slist)
      S = merge_dictionaries!( S, Slist[i] )
      println("md: S: ",S)
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
    ig = chromosome_to_int(g)
    mlist = filter( x->output_values(x)[1]==phenotype, mutate_all( g, funcs, output_chromosomes=true, output_outputs=false ) )
    ihlist = map(h->chromosome_to_int(h),mlist)
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
    println("msd: S: ",S)
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

#=
# Merges int-to-set dictionaries S1 and S2.  The resulting merged dictionary is S1.  S2 is unchanged.
function merge_dictionaries!( S1::Dict{Int64,Set{Int128}}, S2::Dict{Int64,Set{Int128}})
  println("merge_dictionaries: ")
  println("S1: ",S1)
  println("S2: ",S2)
  new_key = maximum(keys(S1)) + 1
  keys2 = [k for k in keys(S2)] 
  println("keys2: ",keys2)
  for ky2 in keys2
    println("S2[",ky2,"]: ",S1[ky2])
    keys1 = [k for k in keys(S1)] 
    println("keys1: ",keys1,"   keys2: ",keys2)
    for ky1 in keys1
      println("S2[",ky2,"]: ",S2[ky2])
      #println("ky2: ",ky2,"  ky1: ",ky1)
      if length(intersect(S1[ky1],S2[ky2])) > 0
        S1[ky1] = union(S1[ky1],S2[ky2])
        println("un S1[",ky1,"]: ",S1[ky1])
        #delete!(S2,ky2)
      else 
        S1[new_key] = deepcopy(S2[ky2])
        println("nk S1[",new_key,"]: ",S1[new_key])
        #delete!(S2,ky2)
        new_key += 1
      end
    end
  end
  S1
end
=#
