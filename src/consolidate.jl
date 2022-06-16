# Consolidates a dataframe with the columns specified below by averaging columns with the same parameter values
# Revised 10/27/20
# Similar to R aggregate function.
function consolidate_dataframe( in_filename::String, out_filename::String; consolidate::Bool=true )
  new_df = DataFrame()
  df = read_dataframe( in_filename )
  #df.n_combined = map(Float64,df.n_combined)  # Only for neutral walk files, remove otherwise
  for n in names(df)  # These two columns must be Float64 in the output dataframe since they are averages
    if n!= "ntries" && n != "evo_count" 
      new_df[!,n] = typeof(df[1,n])[]
    else
      new_df[!,n] = Float64[]
    end
  end
  println("typeof(new_df.evo_count): ",typeof(new_df.evo_count))
  #println("size(new_df): ",size(new_df))
  #println("names(new_df): ",names(new_df))
  #=
  first_pos = findfirst( "0x", df.goallist[1] )[1]
  last_pos = findnext( "]]", df.goallist[1], first_pos )[1]-1
  goals = [parse(UInt16,df.goallist[i][first_pos:last_pos]) for i = 1:length(df.goallist) ]
  df.goal = goals
  =#
  new_names = vcat(["goal"],names(df)[2:end])
  
  println("new_names: ",new_names)
  #=
  for i = 1:size(new_df)[2]
    print(typeof(new_df[1,i])," ")
  end
  println()
  =#
  
  select!(df,new_names)
  df = sort(df,[:goal])
  if !consolidate
    write_dataframe_with_comments( df, in_filename, out_filename )
    return df
  end
  increment = findfirst(x->x!=df.goal[1],df.goal) - 1
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
    #println("size(new_df): ",new_df)
    #println("length(ndf_row): ",length(ndf_row))
    #println("r:",r," ",ndf_row)
    push!(new_df,ndf_row)
  end
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
   write(outfile,"# increement: $increment \n")
   CSV.write( outfile, new_df, append=true, writeheader=true )
  end
  new_df
end
