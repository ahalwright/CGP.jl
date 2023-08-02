using Base.Threads

# Test whether using a vector of Atomic elements elminates any race conditions

# Set up array of Atomic integers and uses threads with Threads.atomic_add!() to increment them
# Gives consistent results showing that this method solves to data race problem
# only slightly slower than non_atomic.
# Example:  test_atomic( 100000 )
function test_atomic( nreps::Int64 )
  len = 20
  test_array = [ Atomic{Int64}(0) for i= 1:len]
  Threads.@threads for i = 1:nreps
    for j = 1:len
      Threads.atomic_add!( test_array[j], 1 )
    end
  end
  [ test_array[i][] for i = 1:len ]
end

# Set up array of integers and uses threads with normal addition increment them
# Gives inconsistent results due to data races
# Example:  test_non_atomic( 100000 )
function test_non_atomic( nreps::Int64 )
  len = 20
  test_array = [ 0 for i= 1:len]
  Threads.@threads for i = 1:nreps
    for j = 1:len
      test_array[j] += 1 
    end
  end
  [ test_array[i] for i = 1:len ]
end

# Set up array of integers and no threads and normal addition increment them
# Faster than either of the above---using threads did not help this application.
# Example:  test_non_threads( 100000 )
function test_non_threads( nreps::Int64 )
  len = 20
  test_array = [ 0 for i= 1:len ]
  for i = 1:nreps
    for j = 1:len
      test_array[j] += 1 
    end
  end
  [ test_array[i] for i = 1:len ]
end
