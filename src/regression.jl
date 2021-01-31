# Functions to compute linear regression statistics.
using GLM, Distributions

function lin_reg( df::DataFrame, Y::Symbol, X::Symbol )
  println("X: ",X,"  Y: ",Y)
  @assert Y in names(df)
  @assert X in names(df)
  d = TDist(size(df)[1]-2)
  global ddf
  ddf = df
  #lreg = lm(@formula( evolvable_count ~ complexity), ddf )
  cmd="lm(@formula($(String(Y)) ~ $(String(X))), ddf )"
  println("cmd: ",cmd)
  lreg = eval(Meta.parse(cmd))
  slope = coef(lreg)[2]
  tempdf = DataFrame( X=>[0.0])
  intercept = predict(lreg,tempdf)[1]
  stderr = stderror(lreg)[2] 
  t = slope/stderr
  pvalue = slope > 0.0 ? ccdf(d,t) : cdf(d,t)
  qminus = quantile(d,0.025)  
  qplus = quantile(d,0.975)   
  lowerci = qplus*stderr + slope   
  upperci = qminus*stderr + slope 
  numinputs = :evo_count in names(df) ? df.numinputs[1] : 3
  ( Y, X, numinputs, slope, intercept, stderr, t, pvalue, lowerci, upperci )
end

function run_lin_reg( df_field_pair_list::Vector{Tuple{DataFrame,Symbol,Symbol}};
    csvfile::String="" )
  lrdf = DataFrame()
  lrdf.Y = Symbol[]
  lrdf.X = Symbol[]
  lrdf.numinputs = Int64[]
  lrdf.slope = Float64[]
  lrdf.intercept = Float64[]
  lrdf.stderr = Float64[]
  lrdf.t = Float64[]
  lrdf.pvalue = Float64[]
  lrdf.lower_ci = Float64[]
  lrdf.upper_ci = Float64[]
  for p in df_field_pair_list
    lrdf_row = lin_reg( p[1], p[2], p[3] )
    push!(lrdf, lrdf_row)
  end
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(lrdf.numinputs[end]))
      CSV.write( f, lrdf, append=true, writeheader=true )
    end
  end
  lrdf
end
