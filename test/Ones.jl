# Construct contexts for up to 7 inputs

const MyInt  = UInt128

function construct_ones( numinputs::Integer )
  Ones = zeros(MyInt,numinputs)
  Ones[1] = 1
  for i = 2:numinputs
    Ones[i] = Ones[i-1] | (Ones[i-1] << 2^(i-2))
  end
  Ones
end

function construct_contexts( numinputs::Integer )
  Ones = construct_ones( numinputs )
  Contexts = [ zeros( MyInt, i ) for i = 1:numinputs ]
  Contexts[1][1] = 0x2
  for i = 2:numinputs
    Contexts[i][1] = Ones[i] << 2^(i-1) 
    for j = 2:i
      #println("(i,j):",(i,j))
      Contexts[i][j] = Contexts[i-1][j-1] << 2^(i-1) | Contexts[i-1][j-1] 
    end
  end
  Contexts
end 
