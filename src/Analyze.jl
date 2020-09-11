using DataFrames
#using StatsBase
using Statistics
using Distributions
using CSV
export spearman_cor, consoidate_dataframe, read_dataframe, write_dataframe_with_comments, write_dataframe


# Returns a dataframe by reading a CSV file with comments that start with "#"
function read_dataframe( df_filename::String )
  df = CSV.read(df_filename,header=true,comment="#")
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

function consolidate_dataframe( in_filename::String, out_filename::String; consolidate::Bool=true )
  new_df = DataFrame()
  new_df.goal=UInt16[]
  new_df.numinputs=Float64[]
  new_df.numoutputs=Float64[]
  new_df.numints=Float64[]
  new_df.levsback=Float64[]
  new_df.ngoals=Float64[]
  new_df.maxsteps=Float64[]
  new_df.gl_reps=Float64[]
  new_df.steps=Float64[]
  new_df.logsteps=Float64[]
  new_df.avgfit=Float64[]
  new_df.nactive=Float64[]
  new_df.redund=Float64[]
  new_df.complex=Float64[]
  new_df.gb_complex=Float64[]
  new_df.degen=Float64[]
  new_df.sdegen=Float64[]
  new_df.f_mutinf=Float64[]
  new_df.mutrobust=Float64[]
  new_df.evolvability=Float64[]     

  df = read_dataframe( in_filename )
  first_pos = findfirst( "0x", df.goallist[1] )[1]
  last_pos = findnext( "]]", df.goallist[1], first_pos )[1]-1
  goals = [parse(UInt16,df.goallist[i][first_pos:last_pos]) for i = 1:length(df.goallist) ]
  df.goal = goals
  new_names = vcat([:goal],names(df)[2:end-1])
  #println("new_names: ",new_names)
  select!(df,new_names)
  df = sort(df,[:goal])
  if !consolidate
    write_dataframe_with_comments( df, in_filename, out_filename )
    return df
  end
  increment = findfirst(x->x!=df.goal[1],df.goal) - 1
  println("increment: ",increment)
  if Int(floor(size(df)[1]/increment)) != Int(ceil(size(df)[1]/increment))
    error("The number of rows in the dataframe must be a multiple of increment")
  end
  for r = 1:increment:Int(floor(size(df)[1]))
    ndf_row = Any[df.goal[r]]
    isum = zeros(Int64,size(df)[2])
    fsum = zeros(Float64,size(df)[2])
    for row = r:(r+increment-1)
      #println("r: ",r,"  row: ",row)
      for i = 2:size(df)[2]
        if typeof(df[row,i]) <: Int64
          isum[i] += df[row,i]
        elseif typeof(df[row,i]) <: Float64 
          fsum[i] += df[row,i]
        else
          error("Illegal type of dataframe element")
        end
      end
    end
    for i = 2:size(df)[2]
      if typeof(df[r,i]) <: Int64
        #print(" ",names(df)[i]," ",isum[i])
        push!(ndf_row,isum[i]/increment)
      elseif typeof(df[r,i]) <: Float64 
        #print(" ",names(df)[i]," ",fsum[i])
        push!(ndf_row,fsum[i]/increment)
      end
      #println()
    end
    #println("r: ",r," ",ndf_row)
    push!(new_df,ndf_row)
  end
  open( out_filename, "w" ) do outfile
    open( in_filename, "r" ) do infile
      line = readline(infile)
      while line[1] == '#'
        #println("line: ",line)
        if line[1:6] != "# cor("
          write(outfile,string(line,"\n"))
        end
        line = readline(infile)
      end
    end
   write(outfile,"# increement: $increment \n")
   CSV.write( outfile, new_df, append=true, writeheader=true )
  end
  new_df
end

function write_dataframe( df::DataFrame, out_filename::String )
  CSV.write( outfile, df, append=true, writeheader=true )
end

# Writes the datarame df to the file out_filename with comments taken from file in_filename.
# No dataframe is read from in_filename.
function write_dataframe_with_comments( df::DataFrame, in_filename::String, out_filename::String )
  open( out_filename, "w" ) do outfile
    open( in_filename, "r" ) do infile
      line = readline(infile)
      while line[1] == '#'
        #println("line: ",line)
        if line[1:6] != "# cor("
          write(outfile,string(line,"\n"))
        end
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

