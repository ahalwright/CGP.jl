using DataFrames
#using StatsBase
using Statistics
using Distributions
using CSV
using StatsBase
using Printf
export spearman_cor, consolidate_dataframe, read_dataframe, write_dataframe_with_comments, write_dataframe
export goal_lookup_from_df, add_edf_col_from_df, frequencies_from_counts, variance1, variance2

# Returns a dataframe by reading a CSV file with comments that start with "#"
function read_dataframe( df_filename::String )
  df = CSV.read(df_filename,DataFrame,header=true,comment="#")
  df
end  

function scor_df( filename::String )
  df = read_dataframe( filename )
  rnames = vcat( [:logsteps], names(df)[12:20])  # result column names
  cdf = DataFrame()
  cdf.names = String.(rnames)
  for i = 1:length(rnames)
    cdf[!,rnames[i]] = zeros(Float64,length(rnames))
    for j = (i+1):length(rnames)
      cdf[j,rnames[i]] = StatsBase.corspearman( df[!,rnames[i]], df[!,rnames[j]] )
    end
  end
  cdf
end

function t_value( v1::Vector{Float64}, v2::Vector{Float64} )
  r = corspearman( v1, v2 )
  r*sqrt((length(v1)-2)/(1-r^2))
end

# Returns a pair of the Spearman correlation between two columns of df and a p-value for this correlation
function spearman_cor( df::DataFrame, name1::Symbol, name2::Symbol )
  v1 = Vector{Float64}(df[!,name1])
  v2 = Vector{Float64}(df[!,name2])
  r = StatsBase.corspearman( v1, v2 )
  t_value = r*sqrt((length(v1)-2)/(1-r^2))
  # p_value is always positive 
  p_value = r>=0.0 ? Distributions.ccdf( TDist(length(v1)-2), t_value ) : Distributions.cdf( TDist(length(v1)-2), t_value ) 
  (r, p_value)
end

function spearman_cor( v1::AbstractVector, v2::AbstractVector )
  r = StatsBase.corspearman( v1, v2 )
  t_value = r*sqrt((length(v1)-2)/(1-r^2))
  # p_value is always positive
  p_value = r>=0.0 ? Distributions.ccdf( TDist(length(v1)-2), t_value ) : Distributions.cdf( TDist(length(v1)-2), t_value )
  (r, p_value)   
end 

# Consolidates a dataframe with the columns specified below by averaging columns with the same parameter values
# Revised 10/27/20
# Similar to R aggregate function.
function consolidate_dataframe( in_filename::String, out_filename::String )
  df = read_dataframe( in_filename )
  cdf = consolidate_dataframe( df )
  open( out_filename, "w" ) do outfile
    open( in_filename, "r" ) do infile
      line = readline(infile)
      while line[1] == '#'
        #println("line: ",line)
        #if line[1:6] != "# cor("
          write(outfile,string(line,"\n"))
        #end
        line = readline(infile)
      end
    end
   #write(outfile,"# consolidate increement: $increment \n")
   CSV.write( outfile, cdf, append=true, writeheader=true )
  end
end

# Consolidates a dataframe with the columns specified below by averaging rows with the same parameter values over all columns except goals
# Assumes that the first column is the "goal" column which may either be Vectors of MyInts or strings.
# The :evo_count column is not averaged since these are cummulative over the rows corresponding to a goal.
# Thus, for the :evo_count column, the last row corresponding to the goal is the consolidated value.
# Revised 10/18/21   Previous values for evo_count were not correct
# Similar to R aggregate function.
function consolidate_dataframe( df::DataFrame )
  new_df = DataFrame()
  #df.n_combined = map(Float64,df.n_combined)  # Only for neutral walk files, remove otherwise
  for n in names(df)
    if n!= "ntries" && n != "evo_count" 
      new_df[!,n] = typeof(df[1,n])[]
    else
      new_df[!,n] = Float64[]
    end
  end
  #println("size(new_df): ",size(new_df))
  #println("names(new_df): ",names(new_df))
  new_names = vcat(["goal"],names(df)[2:end])
  #println("new_names: ",new_names)
  select!(df,new_names)
  df = sort(df,[:goal])
  ff = findfirst(x->x!=df.goal[1],df.goal) 
  increment = ff != nothing ? ff - 1 : size(df)[1]
  #println("increment: ",increment)
  if Int(floor(size(df)[1]/increment)) != Int(ceil(size(df)[1]/increment))
    error("The number of rows in the dataframe must be a multiple of increment")
  end
  for r = 1:increment:Int(floor(size(df)[1]))
    ndf_row = Any[df.goal[r]]   # first entry of row is the current goal
    isum = zeros(Int64,size(df)[2])
    fsum = zeros(Float64,size(df)[2])
    for row = r:(r+increment-1)
      #println("r: ",r,"  row: ",row)
      for i = 2:size(df)[2]
        if names(df)[i] == "evo_count"
          #println("evo_count skip value: ",df[row,:evo_count])
          continue
        elseif typeof(df[row,i]) <: Int64
          isum[i] += df[row,i]
        elseif typeof(df[row,i]) <: Float64 
          fsum[i] += df[row,i]
        else
          error("Illegal type of dataframe element")
        end
      end
    end
    for i = 2:size(df)[2]
      if names(df)[i] == "evo_count"
        #println("evo_count push value: ",df[r+increment-1,:evo_count])
        push!(ndf_row,df[r+increment-1,:evo_count])  # last row for goal has correct evo_count
      elseif typeof(df[r,i]) <: Int64
        #print(" ",names(df)[i]," ",isum[i])
        push!(ndf_row,isum[i]/increment)
      elseif typeof(df[r,i]) <: Float64 
        #print(" ",names(df)[i]," ",fsum[i])
        push!(ndf_row,fsum[i]/increment)
      end
      #println()
    end
    #println("size(new_df): ",new_df)
    #println("length(ndf_row): ",length(ndf_row))
    #println("r:",r," ",ndf_row)
    push!(new_df,ndf_row)
  end
  new_df
end

function write_dataframe( df::DataFrame, out_filename::String )
  CSV.write( out_filename, df, append=true, writeheader=true )
end

# Writes the datarame df to the file out_filename with comments taken from file comments_filename.
# No dataframe is read from comments_filename.
function write_dataframe_with_comments( df::DataFrame, comments_filename::String, out_filename::String )
  open( out_filename, "w" ) do outfile
    open( comments_filename, "r" ) do infile
      line = readline(infile)
      while line[1] == '#'
        #println("line: ",line)
        #if line[1:6] != "# cor("
        #  write(outfile,string(line,"\n"))
        #end
        line = readline(infile)
      end
    end
   #write(outfile,"# increement: $increment \n")
   CSV.write( outfile, df, append=true, writeheader=true )
  end
end

function write_comments( in_filename::String, out_filename::String )
  open( out_filename, "w" ) do outfile
    open( in_filename, "r" ) do infile
      line = readline(infile)
      while line[1] == '#'
        println("line: ",line)
        write(outfile,string(line,"\n"))
        line = readline(infile)
      end
    end
  end
end

# Given a goal represented as a string, extract the df values of the columns fields for that goal.
# Example:  goal_lookup_from_df(df,"0x4208",[:counts10ints]) returns 1192.
function goal_lookup_from_df( df::DataFrame, goal::String, fields::Vector{Symbol} )
  [df[findfirst(x->x==goal,df.goal),f] for f in fields]
end

# Adds a column to dataframe edf whose values are obtained by goal_lookup_from_df( df, goal_string, [field] )
function add_edf_col_from_df( df::DataFrame, edf::DataFrame, field::Symbol )
  edf[!,field] =  zeros(Float64,size(edf)[1])
  for i = 1:size(edf)[1]
    # Extract MyInt goal from edf.goal[i] by parsing.  See Julia documentation for Metaprogramming.
    goal = Meta.parse(edf[i,:goal]).args[2]  
    goal_string = @sprintf("0x%x",goal)   # Convert to string
    edf[i,field] = goal_lookup_from_df(df,goal_string,[field])[1]
  end
  edf
end

# Return the frequencies of the unique elements of V
function frequencies_from_counts( V::Vector{Int64} )
  sort!(V)
  #println("V: ",V) 
  freq = Int64[]
  i = 1
  while i <= length(V)
    push!(freq,1)
    while i < length(V) && V[i+1] == V[i]
      freq[end] += 1
      i += 1
    end
    i += 1
  end
  #println("sum(freq): ",sum(freq),"  freq: ",freq)
  # Two assertions to check correctness
  @assert sum(freq) == length(V)
  @assert length(freq) == length(unique(V))
  freq
end

function pad_missing( V::Vector, new_length::Int64 )
  result = Union{Missing,Int64}[]
    for v in V push!(result,v)
  end
  for i = length(V)+1:new_length
    push!(result,missing)
  end
  result
end
    
# Simple sample variance
function variance1( X::Vector )
  N = length(X)
  avg = sum(X)/N
  sum( (x - avg)^2 for x in X )/(N-1) 
end

# Simple sample variance
function variance2( X::Vector )
  N = length(X)
  ( sum( x^2 for x in X) - sum(X)^2/N )/(N-1)
end
