# Contains functions that are useful in many contexts
export evolvability, run_evolvability, evo_result, test_evo, evo_result_type, run_evolve_g_pairs 
export lg10, set_to_list, write_df_to_csv, MyInt_to_string, string_to_MyInt
export string_to_expression, binomial_prob, binomial_p_value

lg10(x) = iszero(x) ? 0.0 : log10(x)

# Converts a MyInt unsigned integer to a string
# Example:  MyInt_to_string( MyInt(1) ) returns "0x00000001" if MyInt==UInt32
function MyInt_to_string( x::MyInt )
  if MyInt == UInt8
    @sprintf("0x%02x",x)
  elseif MyInt == UInt16
    @sprintf("0x%04x",x)
  elseif MyInt == UInt32
    @sprintf("0x%08x",x)
  elseif MyInt == UInt64
    @sprintf("0x%016x",x)
  elseif MyInt == UInt128
    @sprintf("0x%032x",x)
  else
    error("Illegal type ",MyInt,"  for MyInt in function MyInt_to_string()")
  end
end

# Converts a string representing an MyInt or a string representing a Vector containing a MyInt to a MyInt.
# Example:  string_to_MyInt("0x23")  returns the unsigned integer  0x0023 if MyInt==UInt16
# Example:  string_to_MyInt("UInt16[0x0023]") returns the unsigned integer 0x0023 if MyInt==UInt16
# Example:  string_to_MyInt( "0x004" )  returns 0x00000004 if MyInt==UInt32
function string_to_MyInt( x::AbstractString )
  value = eval(Meta.parse(x))
  if typeof(value) <: Number
    MyInt(value)
  elseif typeof(value) <: Array
    MyInt(value[1])
  else
    error("Illegal argument ",x," to function string_to_MyInt()" )
  end
end

# Converts a string representing a Julia expression to that expression
function string_to_expression( x::AbstractString ) 
  eval(Meta.parse(x))
end

function set_to_list( s::Set )
  [ i for i in s ]
end

# Writes a dataframe to a CSV file
function write_df_to_csv( df::DataFrame, p::Parameters, funcs::Vector{Func}, csvfile::String;
                mutrate::Float64=-1.0, ngens::Int64=-1, popsize::Int64=-1, nreps::Int64=-1, max_steps::Int64=-1,
                numcircuits::Int64=-1, max_mutates::Int64=-1, nsamples::Int64=-1, goal_list::Vector{Vector{MyInt}}=Vector{MyInt}[] )
  open( csvfile, "w" ) do f
    hostname = readchomp(`hostname`)
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )     
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ",funcs)
    if mutrate >= 0.0
      println(f,"# mutrate: ",mutrate)
    end 
    if ngens >= 0
      println(f,"# ngens: ",ngens)
    end 
    if max_steps >= 0
      println(f,"# max_steps: ",max_steps)
    end 
    if nreps >= 0
      println(f,"# nreps: ",@sprintf("%.1e",nreps))
    end 
    if popsize >= 0
      println(f,"# popsize: ",popsize)
    end 
    if numcircuits >= 0
      println(f,"# numcircuits: ",numcircuits)
    end 
    if nsamples >= 0
      println(f,"# nsamples: ",nsamples)
    end 
    if max_mutates >= 0
      println(f,"# max_mutates: ",max_mutates)
    end 
    if length(goal_list) > 0
      println(f,"# goal_list: ",goal_list)
    end 
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end

function assign_to_dataframe_row!( df::DataFrame, index::Int64, row::Vector )
  @assert size(df)[2] == length(row)
  for j = 1:size(df)[2]
    df[index,j] = row[j]
  end
end

# The probabilities argument is a vector of non-negative real numbers which can be normalized into vector of probabilities
# Return an index of this list selected according to these probabilities
# Example: select_from_probabilities([3.,2.,1.,4.]) returns 1 with probability 0.3, 2 with prob 0.2, 3 with prob 0.1, etc.
function select_from_probabilities( probabilities::Vector{<:Number} )
  prob = rand()   # A random number between 0 and 1
  #println("prob: ",prob)
  probabilities = probabilities/sum(probabilities)
  cum_sum = 0.0
  for i = 1:length(probabilities)
    if probabilities[i] < 0
      error("Error in select_from_probabilities:  elements of probabilities must be nonnegative")
    end
    cum_sum += probabilities[i]
    probabilities[i] = cum_sum
  end
  #println("prob: ",prob,"  probabilities: ",probabilities)
  for i = 1:length(probabilities)
    #println("i: ",i,"  probabilities[i]: ", probabilities[i] )
    if prob <= probabilities[i]
      return i
    end
  end
end 

function pheno_filter( probabilities::Vector{<:Number}, p::Parameters, cutoff::Float64, target::Goal )
  # ph_select is a Boolean array over all phenotypes where ph_select[ph[1]] is true if the Hamming dist from ph to target is <= cutoff
  ph_select = map( ph->(hamming_distance( ph, target[1], p.numinputs) <= cutoff), MyInt(0):MyInt(2^2^p.numinputs-1)) 
  new_probabilities = [ (ph_select[i+1] ? probabilities[i+1] : 0.0) for i = MyInt(0):MyInt(2^2^p.numinputs-1) ]
  #for i = 1:length(probabilities)
  #  probabilities[i] = (ph_select[i] ? probabilities[i] : 0.0)
  #end
end

function pheno_filter!( probabilities::Vector{<:Number}, p::Parameters, cutoff::Float64, target::Goal )
  # ph_select is a Boolean array over all phenotypes where ph_select[ph[1]] is true if the Hamming dist from ph to target is <= cutoff
  ph_select = map( ph->(hamming_distance( ph, target[1], p.numinputs) <= cutoff), MyInt(0):MyInt(2^2^p.numinputs-1)) 
  #new_probabilities = [ (ph_select[i+1] ? probabilities[i+1] : 0.0) for i = MyInt(0):MyInt(2^2^p.numinputs-1) ]
  for i = 1:length(probabilities)
    probabilities[i] = (ph_select[i] ? probabilities[i] : 0.0)
  end
end

function test_ph_select( p::Parameters, target::Goal, cutoff::Float64 )
  ph_select = map( ph->(hamming_distance( ph, target[1], p.numinputs) <= cutoff), MyInt(0):MyInt(2^2^p.numinputs-1))
  for i = MyInt(0):MyInt(2^2^p.numinputs-1)
    if ph_select[i+1]
      println( "i: ",i,"  hd: ",hamming_distance(i,target[1],p.numinputs))
    end
  end
end

function pheno_evolve( p::Parameters, funcs::Vector{Func}, start_ph::Goal, target_ph::Goal, pheno_matrix::Matrix{<:Number}; max_reps::Integer=40 )
  current_distance = hamming_distance( start_ph[1], target_ph[1], p.numinputs )
  #println("current_distance: ",current_distance)
  history = [ (start_ph,current_distance) ]
  count = 0
  while count < max_reps &&  current_distance > 0
    probabilities = pheno_matrix[ start_ph[1]+1, : ]
    pheno_filter!( probabilities, p, current_distance, target_ph )
    #println( findall( x->x>0, probabilities ))
    k = select_from_probabilities( probabilities )
    #println("k: ",k,"  probabilities[k]: ",probabilities[k])
    start_ph = [ MyInt(k-1) ]
    current_distance = hamming_distance( start_ph[1], target_ph[1], p.numinputs )
    #println("start_ph: ",start_ph,"  current_distance: ",current_distance)
    push!(history,(start_ph,current_distance))
    count += 1
  end
  history
end
  
# binomial distribution probability for n trials
# Example:  binomial_prob( 10, 0.125, 4, 10 ):  0.02746405079960823 which agrees with https://stattrek.com/online-calculator/binomial
function binomial_prob( n::Int64, p::Float64, lower_bound::Int64, upper_bound::Int64 )
  @assert 0 <= lower_bound <= upper_bound <= n
  sum( k->binomial(n,k) * p^k * (1-p)^(n-k), lower_bound:upper_bound )
end

# One-sided p-value for the one-sided hypothesis that prob <= prob0 where n is the number of trials
# See https://en.wikipedia.org/wiki/Binomial_test
function binomial_p_value( n::Int64, p::Float64, prob0::Float64 )
  k = Int(floor( n*prob0 ))
  binomial_prob( n, p, 0, k )
end
