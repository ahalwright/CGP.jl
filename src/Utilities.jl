# Contains functions that are useful in many contexts
export write_df_to_csv, MyInt_to_string, string_to_MyInt

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
