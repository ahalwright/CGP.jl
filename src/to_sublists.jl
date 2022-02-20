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
