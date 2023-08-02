# Test random_walk---why does not run multiple threads for 4x1
using Base.Threads

# Example: @time test_random_walk(p,28,2^16,100_000_000)   # 704 seconds
function test_random_walk( nthreads::Int64, nphenos::Int64, nreps::Int64 )
  #nphenos = 2^(2^p.numinputs)
  goal_edge_matrix = Array{Atomic{Int64},2}( undef, nphenos, nphenos )
  println("goal_edge_matrix undef allocated")
  Threads.@threads for i = 1:nphenos 
    for j=1:nphenos 
      goal_edge_matrix[i,j]= Atomic{Int64}(0) 
    end 
  end
  println("goal_edge_matrix Atomic initialized")
  Threads.@threads for i = 1:nthreads
    for u = 1:nreps
      j = rand(1:nphenos)
      k = rand(1:nphenos)
      Threads.atomic_add!( goal_edge_matrix[ j, k ], 1 )
    end
  end
  println("goal_edge_matrix updated")
  gem = map(x->x[],goal_edge_matrix)
  println("goal_edge_matrix converted to ints")
  df = matrix_to_dataframe0( gem, collect(MyInt(0):MyInt(nphenos-1)) )
end
  
