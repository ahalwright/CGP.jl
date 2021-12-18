# p = Parameters(3,1,4,5)  # Example
#  funcs = default_funcs(p.numinputs)
#  ec2 = enumerate_circuits( p, funcs); length(ec2) 
#  @time S=find_neutral_components(ec2,0x005a); print_lengths(S)
function find_neutral_components( ch_list::Vector{Chromosome}, phenotype::MyInt )
  p = ch_list[1].params
  funcs = default_funcs(p.numinputs)
  ch_list = filter( x->output_values(x)[1]==phenotype, ch_list )
  if length(ch_list) == 0
    error("no genotypes that map to the given phenothype ",[phenotype])
  end
  println("length(ch_list): ",length(ch_list))
  readline()
  S = Dict{Int64,Set{Int128}}()
  new_key = 1
  for g in ch_list
    ig = chromosome_to_int(g)
    mlist = filter( x->output_values(x)[1]==phenotype, mutate_all( g, funcs, output_chromosomes=true, output_outputs=false ) )
    #println("mlist[1]: ",mlist[1])
    ihlist = map(h->chromosome_to_int(h),mlist)
    push!(ihlist,ig)
    ihset = Set(ihlist)
    print("ig: ",ig,"  lenght(ihset): ",length(ihset),"   ")
    if length(ihset) > 0
      for ky in keys(S)
        #println("ky: ",ky,"  S[ky]: ",S[ky])
        if length( intersect( ihset, S[ky] ) ) > 0
          union!( ihset, S[ky] )
          delete!( S, ky )
        end
      end
      S[ new_key ] = ihset
      println("length(S[",new_key,"]) = ",length(S[new_key]))
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

function print_lengths(S)
  for ky in keys(S) 
    println("ky: ",ky,"  length(S[ky])): ",length(S[ky])) 
  end
end
#=
function chromosome_to_int( ch::Chromosome, funcs::Vector{Func}=default_funcs(ch.params.numinputs); maxarity::Int64=2 )
  result = Int128(0)
  for i = 1:length(ch.interiors)
    #println("i: ",i)
    result += gate_int( i, ch, maxarity, funcs )
    multiplier = i < length(ch.interiors) ? length(funcs)*(ch.params.numinputs+i+1)^maxarity : 1
    result *= multiplier
    #result = i < length(ch.interiors) ? result*length(funcs)*(ch.params.numinputs+i+1)^maxarity : result
    #println("i: ",i,"  gate_int: ",gate_int( i, ch, maxarity, funcs ),"  multiplier: ",multiplier, "  result: ",result)
    #println("ci i: ",i,"  gate_int: ",gate_int( i, ch, maxarity, funcs ),"  result: ",result)
  end
  #println("result: ",result)
  result
end

function gate_int( i::Int64, ch::Chromosome, maxarity::Int64, funcs::Vector{Func} )
  #println("i: ",i)
  #println("ch.interiors[i].func: ",ch.interiors[i].func)
  func_list = [ f.func for f in funcs ]   # necessary because the == operator doesn't work on structs
  numinputs = ch.params.numinputs
  funcs_int = findfirst(x->x==ch.interiors[i].func.func,func_list)[1]-1
  gate_inputs = ch.interiors[i].inputs
  il = inputsList( maxarity, numinputs+i )
  ni = length(il)   # multiplier which should be (numinputs+i)^maxarity
 #println("ni: ",ni,"  maxarity: ",maxarity,"  numinputs+i: ",numinputs+i)
  @assert ni == (numinputs+i)^maxarity
  inputs_int = findfirst(x->x==gate_inputs,il)
  funcs_int*ni + inputs_int
end

function inputsList( numinputs::Int64, maxval::Int64 )
  if numinputs == 1
    return [ [i] for i = 1:maxval ]
  end
  result = inputsList( numinputs-1, maxval )
  new_result = Vector{Int64}[]
  for r in result
    for i = 1:maxval
      dcr = deepcopy(r)
      push!(dcr,i)
      push!(new_result,dcr)
      #println("i: ",i,"  dcr: ",dcr)
    end
  end
  #println("numimnputs: ",numinputs,"  new_result: ",new_result)
  new_result
end      
=#
