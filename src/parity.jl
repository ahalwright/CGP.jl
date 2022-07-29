# even_parity() computes the goal for evolution of a circuit to solve the even parity problem
# odd_parity() computes the goal for evolution of a circuit to solve the odd parity problem
function odd_parity(numinputs::Int64)
  gb = reverse(get_bits( construct_context(numinputs),numinputs))
  result = MyInt(0)
  for gbi in gb
    result <<= 1
    #result |= MyInt(1 - (count_ones(gbi) % 2))
    result |= MyInt(1 - (count_ones(gbi) & 1))
  end
  result
end

function even_parity(numinputs::Int64)
  gb = reverse(get_bits( construct_context(numinputs),numinputs))
  result = MyInt(0)
  for gbi in gb
    result <<= 1
    result |= MyInt(count_ones(gbi) & 1)
  end
  result
end

even(x::Integer) = Int64(1-(x & MyInt(1))) 
odd(x::Integer) = Int64(x & MyInt(1))

function test_odd_parity( numinputs )
  blist = [  MyInt_to_binary_list(MyInt(x),numinputs) for x =  MyInt(2^numinputs-1):-1:0x0 ]
  clist = map(x->even(sum(x)),blist)
  to_bin = to_binary(even_parity(numinputs),2^numinputs)   
  @assert to_bin == clist   
end  

function test_even_parity( numinputs )
  blist = [  MyInt_to_binary_list(MyInt(x),numinputs) for x =  MyInt(2^numinputs-1):-1:0x0 ]
  clist = map(x->odd(sum(x)),blist)
  to_bin = to_binary(odd_parity(numinputs),2^numinputs)   
  @assert to_bin == clist   
end  

function MyInt_to_binary_list( x::MyInt, size::Int64 )
  #size = 64  # If the type of x is UInt64, must be 64.  If the type of x is UInt32, must be 32.
  if MyInt==UInt16
    @assert size <= 4
  end
  if MyInt==UInt32 
    @assert size <= 5
  end
  if MyInt==UInt64 
    @assert size <= 6
  end
  if MyInt==UInt128
    @assert size <= 7
  end
  #println("size: ",size)
  shift = size-1
  result = fill(MyInt(0),size)
  for i = 1:size
    result[i] = (x >> shift) & convert(MyInt,1)
    shift -= 1
  end
  result
end

function shift_extend( p::Parameters, maxsteps::Int64; parity::Bool=:true, extend::Int64=3 )
  funcs = default_funcs(p)
  c = random_chromosome(p)
  if parity
    goal = [even_parity(p.numinputs)]
  else
    goal= [rand_shift_extend(extend,p.numinputs)]
  end
  #println("goal: ",goal)
  (nc,steps)=neutral_evolution(c,funcs,goal,maxsteps)
  ssteps =  steps < maxsteps ? steps : 0
  failure = ssteps == 0 ? 1 : 0
  (parity,p.numinputs,p.numinteriors,p.numlevelsback,ssteps,failure)
end

function run_shift_extend( numinputs::IntRange, reps::Int64, maxsteps::Int64; 
    csvfile::String="", parity::Bool=:true, extend::Int64=3 )
  number_gates_multiplier = 5
  #number_gates_multiplier = 7  # Changed from 5 to 7 on 4/29 to do experiment with 15 inputs
  df = DataFrame()
  df.parity = Bool[]
  df.numinputs = Int64[]
  df.numgates = Int64[]
  df.levelsback = Int64[]
  df.repetitions = Int64[]
  if !parity
    df.extend = Int64[]
  end
  df.average_steps = Float64[]
  df.median_steps = Float64[]
  for ni in numinputs
    if extend >= ni
      error("The keywork parameter extend=",extend,"  must be less than numinputsi=",ni)
    end  
    ng = number_gates_multiplier*ni
    levelsback = Int(floor(ng/2))
    p = Parameters( ni, 1, ng, levelsback )
    funcs = default_funcs(p.numinputs)
    rows = pmap(x->shift_extend(p,maxsteps,parity=parity,extend=extend),collect(1:reps))
    sum_steps = 0
    sum_failures = 0
    steps_list = Int64[]
    for r in rows
      if r[5] > 0
        push!(steps_list,r[5])
      else
        sum_failures += r[6]
      end
    end
    r = rows[1]
    if parity
      row = (r[1],r[2],r[3],r[4],reps-sum_failures,mean(steps_list),median(steps_list))
    else
      row = (r[1],r[2],r[3],r[4],reps-sum_failures,extend,mean(steps_list),median(steps_list))
    end 
    push!(df,row)
  end
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(numinputs[end]))
      println(f,"# MyInt: ", MyInt )
      println(f,"# reps: ",reps)
      println(f,"# extend: ",extend)
      println(f,"# maxsteps: ", maxsteps )
      println(f,"# number_gates_multiplier: ",number_gates_multiplier)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

function reverse_bits( x::MyInt, size::Int64 )
  result = MyInt(0)
  mask = MyInt(1)
  final = MyInt(1) << (size-1)
  result |= mask & x != MyInt(0) ? MyInt(1) : MyInt(0)
  while mask != final
    result <<= 1
    mask <<= 1
    result |= mask & x != MyInt(0) ? MyInt(1) : MyInt(0)
    #@printf("mask: 0x%x  result: 0x%x\n",mask,result)
  end
  result
end

function rand_shift_extend( orig_nbits::Int64, result_nbits::Int64; parity::Bool=:false )
  @assert result_nbits > orig_nbits
  if parity 
    sgl = MyInt(0x96)
  else
    sgl = randgoal(orig_nbits,1)[1]
  end
  #@printf("sgl: 0x%x\n",sgl)
  rgl = sgl |= sgl | (reverse_bits(sgl,2^orig_nbits) << 2^orig_nbits)
  #@printf("rgl: 0x%x\n",rgl)
  for i = (orig_nbits+2):result_nbits
    rgl |= (rgl << 2^(i-1))
    #@printf("i: %d  rgl: 0x%x\n",i,rgl)
  end
  rgl
end
