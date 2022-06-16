# Produces a list of all circuits corresponding the paramters p 
function enumeratecircuits_ch( p::Parameters, funcs::Vector{Func}=default_funcs(p.numinputs); maxarity::Int64=2 )
  enumeratecircuits_ch( p, p.numinteriors, funcs, maxarity=maxarity )
end

# Recursive helper function
function enumeratecircuits_ch( p::Parameters, numints::Int64, funcs::Vector{Func}; maxarity::Int64=2 )
  #println("enum_circs_rec: numints: ",numints)
  if numints == 0
    result = [Chromosome( p, map(i->InputNode(i),collect(1:p.numinputs)), InteriorNode[], [OutputNode(p.numinputs)], 0.0, 0.0 )]
    return result
  end
  prev_result = enumeratecircuits_ch( p, numints-1, funcs ) 
  result = Chromosome[]
  for prev_ch in prev_result
    #println(prev_ch)
    for inputs in inputs_List( maxarity, max(1,p.numinputs+numints-p.numlevelsback), p.numinputs+numints-1 )
      #println("p.numinputs+numints-1: ",p.numinputs+numints-1)
      for func in funcs
        new_interiors = vcat(prev_ch.interiors,[InteriorNode(func,inputs)])
        new_ch = Chromosome( p, prev_ch.inputs,new_interiors,[OutputNode(p.numinputs+numints)],0.0,0.0) 
        push!(result,new_ch)
      end
    end
  end
  result
end

# In all non-recursive calls to this function, numinputs is maxarity
function inputs_List( numinputs::Int64, minval::Int64, maxval::Int64 )
  #println("inputsList numimnputs: ",numinputs,"  minval: ",minval,"  maxval: ",maxval )
  if numinputs == 1
    return [ [i] for i = minval:maxval ]
  end
  result = inputs_List( numinputs-1, minval, maxval )
  new_result = Vector{Int64}[]
  for r in result
    for i = minval:maxval
      dcr = deepcopy(r)
      push!(dcr,i)
      push!(new_result,dcr)
      #println("i: ",i,"  dcr: ",dcr)
    end
  end
  #println("numimnputs: ",numinputs,"  input_List_result: ",new_result)
  new_result
end     

