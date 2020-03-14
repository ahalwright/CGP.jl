# Construct contexts for up to 7 inputs if MyInt = UInt128
# Example of context for 2 inputs
# Context[1] = 0x0c = 1100
# Context[2] = 0x0a = 1010
# We consider this as a 2 by 4 bit matrix where Context[i] is the ith row.
# Each column of this matrx corresponds to 2 inputs to the circuit.
# All columns of the matrix correspond to all possible inputs to the circuit.
# Column 1 corresponds to input1 = 1, input2 = 1
# Column 2 corresponds to input1 = 1, input2 = 0
# Column 3 corresponds to input1 = 0, input2 = 1
# Column 4 corresponds to input1 = 0, input2 = 0
# Thus, the 4-bit output string describes the output of the circuit for all inputs.
# This specifies a truth table where columns of the Context matrix correspoing to
#   rows of the truth table.

# When there are 3 inputs, Context has 3 components, one for each input. 
# Thus, the matrix has 3 rows and 8 = 2^3 columns and the function computed by the
#   circuit is specified by the 8 bit output string.

# The logic function multi-output circuit is specified by sequence of output bit 
#   strings for each output, or by the concatenation of these bit strings into 
#   a single bit string.

export construct_ones, construct_contexts

#const MyInt  = UInt8

# Ones[i] is a "mask" with 1 bits in positions 1 to 2^i where position 1 is least significant
function construct_ones( numinputs::Integer )
  Ones = zeros(MyInt,numinputs)
  Ones[1] = 3
  for i = 2:numinputs
    Ones[i] = Ones[i-1] | (Ones[i-1] << 2^(i-1))
  end
  Ones
end

function construct_contexts( numinputs::Integer )
  Ones = construct_ones( numinputs )
  Contexts = [ zeros( MyInt, i ) for i = 1:numinputs ]
  Contexts[1][1] = 0x2
  for i = 2:numinputs
    Contexts[i][1] = Ones[i-1] << 2^(i-1) 
    for j = 2:i
      #println("(i,j):",(i,j))
      Contexts[i][j] = Contexts[i-1][j-1] << 2^(i-1) | Contexts[i-1][j-1] 
    end
  end
  Contexts
end 
