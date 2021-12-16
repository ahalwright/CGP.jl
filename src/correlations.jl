using DataFrames, CSV, StatsBase
export correlations, scatter_plot
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
  println("delete_inds: ",delete_inds)
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
  #=
  for i = 1:sz-1
    for j = i+1:sz
      sc = spearman_cor( df[i,:values], df[j,:values] )
      #println("(i,j): ",(i,j),"  sc[1]: ",sc[1],"  sc[2]: ",sc[2])
      cor_matrix[i,j] = sc[1]
      pval_matrix[i,j] = sc[2]
    end
  end
  =#
  cor_df = matrix_to_dataframe( cor_matrix, df )
  pval_df = matrix_to_dataframe( pval_matrix, df )
  #(df,cor_df,pval_df)
  df
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
 
# cshort_name is the variable that determimes marker color and size
function scatter_plot( adf::DataFrame, xshort_name::String, yshort_name::String, cshort_name::String=""; data_subdir::String="10_29_21" )
  if !(xshort_name in adf.short_name)
    error("xshort_name ",xshort_name," is not in the dataframe in scatter_plot()")
  end
  if !(yshort_name in adf.short_name)
    error("yshort_name ",yshort_name," is not in the dataframe in scatter_plot()")
  end
  if length(cshort_name)>0 && !(cshort_name in adf.short_name)
    error("cshort_name ",cshort_name," is not in the dataframe in scatter_plot()")
  end
  if length(cshort_name)>0 
    cvalues = adf[adf.short_name.==cshort_name,:values][1]
    cmin=findmin(cvalues)[1]
    cmax=findmax(cvalues)[1]
    # Normalized values
    nvalues = map(x->(1.0/(cmax-cmin))*(x-cmin),cvalues)
    color_values = map(x->RGB(x,0.0,1.0-x),nvalues)  # Large red shows high values, small blue shows low values
    size_values = 5 .+ (10 .* nvalues)
    scatter( adf[ adf.short_name.==xshort_name,:values], adf[ adf.short_name.==yshort_name,:values], smooth=true, c=color_values, markersize=size_values )
  else
    scatter( adf[ adf.short_name.==xshort_name,:values], adf[ adf.short_name.==yshort_name,:values], smooth=true )
  end
  scatter!( xlabel=xshort_name, ylabel=yshort_name )
  xparams = xshort_name[1:2] == "CR" ? "8_5" : "8_2"
  yparams = yshort_name[1:2] == "CR" ? "8_5" : "8_2"
  params = xparams==yparams ? xparams : yparams * xparams
  if length(cshort_name) == 0
    title="$yshort_name vs $xshort_name 3x1 $params"
  else
    title="$yshort_name vs $xshort_name markers $cshort_name 3x1 $params"
  end
  scatter!( title=title, legend=:none )
  #readline()
  # current working directory (pwd()) must be the data subdirectory
  fname = data_subdir * "/" * title * ".png"
  println("fname: ",fname)
  savefig( fname )
end
