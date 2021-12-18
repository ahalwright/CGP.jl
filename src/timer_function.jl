using CSV
using Dates

function timer( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, csvfile )
  println(default_funcs(2))
  result =
    @timed run_geno_pheno_evolution( iterations, numinputs, numoutputs, numinteriors, ngoals, maxsteps, levelsback, csvfile )
  df = result[1]
  ttime = result[2]
  println("run time in minutes: ",ttime/60)
  hostname = chomp(open("/etc/hostname") do f read(f,String) end)
  open( csvfile, "w" ) do f
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    println(f,"# run time in minutes: ",ttime/60)
    println(f,"# funcs: ", Main.CGP.default_funcs(numinputs[end]))
    #println(f,"# nodearity: ",nodearity)
    #println(f,"# active_only: ",active_only)
    CSV.write( f, df, append=true, writeheader=true )
  end
  df
end


