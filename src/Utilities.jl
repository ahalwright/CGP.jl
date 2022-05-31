export write_df_to_csv
# Writes a dataframe to a CSV file
function write_df_to_csv( df::DataFrame, p::Parameters, funcs::Vector{Func}, csvfile::String;
    mutrate::Float64=-1.0, ngens::Int64=-1, popsize::Int64=-1, nreps::Int64=-1, max_steps::Int64=-1,
    numcircuits::Int64=-1, max_mutates::Int64=-1, goal_list::Vector{Vector{MyInt}}=Vector{MyInt}[] )
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end) 
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
