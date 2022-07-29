function find_desirable_goals( p::Parameters, goal::MyInt, mindist::Int64, tries::Int64 )
  desirable_list = MyInt[]
  count = 0
  for i = 1:tries
    if 2^p.numinputs < MyIntBits(MyInt)
      g = rand(MyInt(0):MyInt(2^2^p.numinputs-1))
    else
      g = rand(MyInt(0):typemax(MyInt))
    end
    if hamming( g, goal ) <= mindist
      count += 1
      push!(desirable_list,g)
    end
  end
  (desirable_list,count)
end

# Find complex goals by screening goals by difficulty of evolution
function find_complex_goals( p::Parameters, funcs::Vector{Func}, nevolutions::Int64, maxsteps::Int64, tries::Int64 )
  goal_list = Tuple{MyInt,Int64}[]
  for i = 1:tries
    if 2^p.numinputs < MyIntBits(MyInt)
      g = rand(MyInt(0):MyInt(2^2^p.numinputs-1))
    else
      g = rand(MyInt(0):typemax(MyInt))
    end
    sum_dist = 0.0
    for j = 1:nevolutions
      c = random_chromosome(p,funcs)
      (nc,steps) = neutral_evolution( c, funcs, [g], maxsteps, funcs=funcs )
      sum_dist += hamming( output_values(nc)[1], g )
    end
    push!(goal_list,(g,sum_dist))
  end
  goal_list
end
