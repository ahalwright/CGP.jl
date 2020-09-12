# Creates an array of the outputs created by nreps random chromosomes.
# Prints (shows) this array, showing elemets greater than a show_increment parameter, or less than or equal to  show_increment.
# Or shows the count of outputs satisfying these conditions.
# Since an array is used to store the outputs, numinputs for chromosomes is limited to 4.
# In order to handle larger numinputs, a hash table would need to be used---not implemented at this time.
# Example:
# julia> ol = count_outputs_parallel( 100, 2, 1, 4, 4 );
# julia> show_output_list( ol, 2 )  # shows outputs created more than 2 times
# output: 0X0:  24
# output: 0X1:  8
# output: 0X2:  3
# output: 0X3:  4
# output: 0X5:  11
# output: 0X6:  3
# output: 0X7:  6
# output: 0X8:  5
# output: 0X9:  3
# output: 0Xa:  13
# output: 0Xb:  7
# output: 0Xc:  6
# output: 0Xd:  4
# output: 0Xe:  8
# output: 0Xf:  18
# count of shown outputs: 15
# julia> show_outputs_list( ol, 2, show_small=true )  # shows outputs created at most twice
# output: 0X4:  2
# count of shown outputs: 1
# See notes/diargy9_9_20.txt for more examples.
using DataFrames
using Distributed
using Printf
export count_outputs, count_outputs_parallel, write_to_file, read_file, show_outputs_list, read_counts_files
export add_counts_to_dataframe

MyFunc = Main.CGP.MyFunc

# Creates the array of output counts
function create_count_outputs_list( numinputs::Integer, numoutputs::Integer )
  return fill( convert(MyFunc,0), numoutputs*2^2^numinputs )
end

# increments the array of output counts
function increment_count_outputs_list( output::Goal, outputs_list::Vector{MyFunc}, numinputs::Int64 )
  outputs_list[concatenate_outputs(output,numinputs)+1] += 1
end

# If numoutputs > 1, combines the outputs into a single unsigned integer of type MyFunc
# If numoutputs == 1, extract the one compponent
function concatenate_outputs( output::Goal, numinputs::Int64 )
  if length(output) == 1
    return MyFunc(output[1])
  else
    result = MyFunc(0)
    for i = 1:length(output)
      result <<= 2^numinputs
      result |= output[i]
    end
    result
  end
end

# Return an output list of the number of times that an output was produced by randomly generating chromosomes with these parameters
function count_outputs( nreps::Int64, numinputs::Integer, numoutputs::Integer, numinteriors::Int64, numlevelsback::Integer )
  p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numinteriors, numlevelsback=numlevelsback ) 
  funcs = default_funcs(numinputs)
  outlist = create_count_outputs_list( numinputs, numoutputs )
  for _ = 1:nreps
    c = random_chromosome( p, funcs )
    output = output_values( c )
    increment_count_outputs_list( output, outlist, numinputs )
  end
  outlist
end

# Return an output list of the number of times that an output was produced by randomly generating chromosomes with these parameters
function count_outputs_parallel( nreps::Int64, numinputs::Integer, numoutputs::Integer, numinteriors::Int64, numlevelsback::Integer ) 
  if nprocs() > 1
    nreps_p = Int(round(nreps/(nprocs()-1)))
  else
    nreps_p = 1
  end
  println("nreps_p: ",nreps_p)
  #mapreduce(x->count_outputs( nreps_p, numinputs, numoutputs, numinteriors, numlevelsback ), +, collect(1:nprocs()))
  reduce(+, pmap( x->count_outputs( nreps_p, numinputs, numoutputs, numinteriors, numlevelsback ), collect(1:nprocs())))
end

# print the elements of the output list that are greater than increment
# if show_small==true, print the elements of the output list that are less than or equal to  increment
function show_outputs_list( outputs_list::Vector{MyFunc}, show_increment::Int64=0; show_small::Bool=false, count_only::Bool=false )
  count = 0
  for i = 0:length(outputs_list)-1
    if !show_small
      if outputs_list[i+1] > show_increment 
        if !count_only
          @printf("output: 0X%x:  %d\n",i,outputs_list[i+1])
        end
        count += 1
      end
    else
      if outputs_list[i+1] <= show_increment 
        if !count_only
          @printf("output: 0X%x:  %d\n",i,outputs_list[i+1])
        end
        count += 1
      end
    end
  end
  println("count of shown elements: ",count)
end

# Write output counts list to a file.  
# Optional list of comments which are written to a file preceded by a # character.
# Example:  write_to_file( ol, "../data/9_8/test.csv", "10^6 reps Raman funcs", "3 inputs, 1 output, 7 interiors, 4 levsback",hex=true )
function write_to_file( outputs_list::Vector{MyFunc}, filename::String, comments::String... ; 
      hex::Bool=false )
  open( filename, "w" ) do f
    for s in comments
      println(f,"# ",s)
    end
    for i = 1:length(outputs_list)
      if hex
        @printf(f,"0x%x\n",outputs_list[i])
      else
        println(f, outputs_list[i] )
      end
    end
  end
end

# Read a counts file written by write_to_file(). 
# Returns a list whose elements are of type out_type.
# Options for out_type include Int64, UInt64, UInt32,
# Lines starting with "#" are comments
function read_file( filename::String, out_type::Type )
  result = out_type[]
  i = 0
  for row in CSV.File(filename, header=false, comment="#" )
    #println("i: ",i,"  ",row.Column1)
    push!(result,parse(out_type,row.Column1))
    i += 1
  end
  result
end

# Read multiple counts files.   
# Returns a 2-tuple where the first element is a list of goals written in hex form,
#    and the second is a vector of counts lists, one for each file read.
function read_counts_files(filename::String...)
  my_goals = Vector{String}[]
  my_counts = Vector{Int64}[]
  i=1
  for fname in filename
    open(fname) do f
      lines = readlines(f)
      if i > 1
        @assert length(lines) == length(my_goals[i-1])
      end
      filter!(x->x[1]!='#',lines)
      goals = [ @sprintf("0x%x",i) for i = 0:(length(lines)-1) ]
      if i > 1
        @assert goals == my_goals[i-1]
      end
      push!(my_goals,goals)
      counts = [ parse(Int64, lines[i+1] ) for i = 0:(length(lines)-1)]
      push!(my_counts,counts)
    end
  end
  (my_goals[1],my_counts)
end

# Add multiple counts vectors as additional columns to dataframe df.
# The hames of the columns are given in the names vector.
# The filename to read the counts vector is filename.
# The dataframe can be read using function read_dataframe() in Analyze.jl
function add_counts_to_dataframe( df::DataFrame, names::Vector{Symbol}, filename::String... )
  (goals,counts) = read_counts_files( filename... )
  df.goal = goals
  i = 1
  for nm in names
    insert!(df, size(df)[2]+1, counts[i], nm )
  end
  df
end

