# test that sampling approximates weight evolvability.

function test_weighted_evolvabilty( nsamples )
  counts = zeros(Int64,8)
  ph = Vector{Int64}[]
  push!(ph,[1,1,2,4,0])
  push!(ph,[1,1,2,3,0])
  push!(ph,[1,4,2,0,0])
  push!(ph,[1,6,2,0,0])
  push!(ph,[1,3,5,0,0])
  println(ph)
  samples = rand(ph,nsamples)
  for s in samples
    for c in s
      if c > 0
        counts[c] += 1
      end
    end
  end
  counts
end

