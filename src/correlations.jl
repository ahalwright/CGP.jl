using DataFrames, CSV, StatsBase
# read a csv which includes vectors of interest, and compute all cross correlations
function correlations( csvfile::String, row_list::Vector{Int64}=Int64[] )
  #df = CSV.read("../notes/correlation_csvs/corr.csv",DataFrame,comment="#")
  df = CSV.read(csvfile,DataFrame,comment="#")
  row_list = length(row_list)==0 ? (1:size(df)[1]) : row_list
  println("row_list: ",row_list)
  sz = length(row_list)
  if Sys.islinux()
    cd("../data")
  elseif Sys.iswindows()
    cd("../../complexity/data")
  end
  value_vects = Vector[]
  for i = row_list
    tmpdf = read_dataframe(df[i,"CSV"]);
    col = tmpdf[:,df[i,"colname"]]
    if df[i,"Dataframe"]=="LCS"   # temporary fix for NaN values
      println("i: ",i,"   name: ",df[i,"Dataframe"],"  colname: ",df[i,"colname"],"  bad value: ",col[108])
      col[108] = 4.0  # Approximate mean complexity
      col[135] = 4.0  # Approximate mean complexity
    end
    #=
    if df[i,"Dataframe"]=="CRE" && df[i,"colname"]=="robustness"  # temporary fix for NaN value
      println("i: ",i,"   name: ",df[i,"Dataframe"],"  colname: ",df[i,"colname"],"  bad value: ",col[151])
      col[151] = 0.2  # a guess
    end
    =#
    push!(value_vects,col)
  end
  delete_inds = setdiff( collect(1:size(df)[1]), row_list )
  delete!(df,delete_inds)
  println("df.short_name: ",df.short_name)
  if Sys.islinux()
    cd("../src")
  elseif Sys.iswindows()
    cd("../../../src")
  end
  insertcols!(df,:values=>value_vects)
  #df
  cor_matrix = fill(0.0,sz,sz)
  pval_matrix = fill(0.0,sz,sz)
  for i = 1:sz-1
    for j = i+1:sz
      sc = spearman_cor( df[i,:values], df[j,:values] )
      #println("(i,j): ",(i,j),"  sc[1]: ",sc[1],"  sc[2]: ",sc[2])
      cor_matrix[i,j] = sc[1]
      pval_matrix[i,j] = sc[2]
    end
  end
  cor_df = matrix_to_dataframe( cor_matrix, df )
  pval_df = matrix_to_dataframe( pval_matrix, df )
  (df,cor_df,pval_df)
end

function matrix_to_dataframe( matrix::Matrix{Float64},df::DataFrame)
  @assert size(matrix)[1] == size(matrix)[2]
  rmatrix = map(x->round(x,sigdigits=3),matrix)
  mdf = DataFrame()
  insertcols!(mdf, :name=>df.short_name[1:(end-1)]) 
  for i = 2:size(rmatrix)[1]
    insertcols!(mdf, Symbol(df.short_name[i]) => rmatrix[1:(end-1),i] )
  end 
  mdf
end
 
