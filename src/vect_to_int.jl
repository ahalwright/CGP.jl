# Convert a vectors of integers to a unique integer.
# The components of the vector are 1-based rather than 0-based.

# Example 1:  
# julia> vect_to_int( [3,4,5,6],[10,10,10,10])
#  2345
# Example 2:
# julia> vect_to_int([3,2,1,4],[7,2,5,8])
#  203

function vect_to_int( ints::Vector{Int64}, max::Vector{Int64} )
  result = Int128(ints[1])-1
  for i = 2:length(ints)
    result *= max[i]
    result += ints[i]-1
  end
  result
end

# julia> int_to_vect( 2345, [10,10,10,10])
# 4-element Array{Int64,1}:
#  Int64[ 3, 4, 5, 6 ]
# julia> int_to_vect(203,[7,2,5,8])
# 4-element Array{Int64,1}:
#  Int64[ 3, 2, 1, 4] 

function int_to_vect( int::Integer, max::Vector{Int64} )
  ints = zeros(Int64,length(max))
  for i = length(max):-1:1
    ints[i] = int % max[i] + 1
    int = div(int,max[i])
  end
  ints
end
    
