
function test_circuit_codes( p::Parameters )
  funcs = default_funcs( p.numinputs )
  c = random_chromosome( p, funcs )
  cc = circuit_code( c )
  # See https://en.wikibooks.org/wiki/Introducing_Julia/Strings_and_characters#Streams
  iobuffer = IOBuffer()   
  print_build_chromosome(iobuffer,c)
  c_str = String(take!(iobuffer))
  ctc = code_to_circuit(cc,p)
  iobuffer = IOBuffer()
  print_build_chromosome(iobuffer,c)
  ctc_str = String(take!(iobuffer))
  if c_str != ctc_str
    println("error c != ctc ")
    println("c_str:   ",c_str)
    println("ctc_str: ",ctc_str)
  end
end

function ttt(j)
  for i = 1:4
   println(j)
   r = j % 2
   j ÷= 5
   println("i: ",i,"  r: ",r,"  j: ",j)
   r = j % 2
   j ÷= 2
   println("i: ",i,"  r: ",r,"  j: ",j)
   r = j % 5
   j ÷= 2
   println("i: ",i,"  r: ",r,"  j: ",j)
  end
end

function intc( j::Integer, p::Parameters )
  result = Int64[]
  for i = 1:p.numinteriors
    jmod = j % 2
    push!( result, jmod )
    println("i: ",i,"  j % 2: ",jmod,"  j: ",j)
    j ÷= 2
    jmod = j % 2
    push!( result, jmod )
    println("i: ",i,"  j % 2: ",jmod,"  j: ",j)
    j ÷= 2
    jmod = j % 5
    push!( result, jmod )
    println("i: ",i,"  j % 5: ",jmod,"  j: ",j)
    j ÷= 5
    println("res: ",result)
  end
  transpose(result)
end 
function int_to_circuit_code( c_int::Integer, p::Parameters )
  c_int = Int128(c_int)
  c_code = zeros(Int64,3*p.numinteriors)
  k = 3*p.numinteriors
  for i = p.numinteriors:-1:1
    multiplier = min(p.numlevelsback,i-1+p.numinputs)
    println("i: ",i,"  multiplier: ",multiplier)
    for j = p.nodearity:-1:1
      c_int_mod= c_int % multiplier
      c_int ÷= multiplier
      println("i: ",i,"  j: ",j,"  c_int: ",c_int,"  c_int_mod: ",c_int_mod)
      c_code[k] = c_int_mod
      k -= 1
    end
    multiplier = length(funcs)
    println("i: ",i,"  multiplier: ",multiplier)
    c_int_mod= c_int % multiplier
    c_int ÷= multiplier
    println("i: ",i,"  c_int: ",c_int,"  c_int_mod: ",c_int_mod)
    c_code[k] = c_int_mod
    k -= 1
  end
  c_code
end          
