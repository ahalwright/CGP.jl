# Returns the number active not inluding inputs.
function num_active_lc( circ::LinCircuit, funcs::Vector{Func}=default_funcs(circ.params) )
  p = circ.params
  ninstructions = p.numinteriors
  res = num_active_lc( circ, numinstructions, numinstructions+1, funcs )
  reduce(+,res[2])
end

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
