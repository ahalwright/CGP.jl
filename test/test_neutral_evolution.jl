# Test how well neutral evolution works with different parameter settings.

function run_test_neutral_evolution( repetitions::Int64, numinputs::Int64, numinstructions_rng::IntRange, numregisters_rng::IntRange, max_steps;
      cartesian::Bool=:false, csvfile::String="" )
  if nprocs() > 2
    process_repetitions = Int(floor(repetitions/(nprocs()-1)))
    println("process_repetitions: ",process_repetitions)
    df_list = pmap( x->test_neutral_evolution( process_repetitions, numinputs, numinstructions_rng, numregisters_rng, max_steps, cartesian=cartesian ), 
        collect(1:(nprocs()-1)) )
    result_df = df_list[1]
    for i = 2:length(df_list)
      result_df = average_compatible_dataframes( result_df, df_list[i] )
    end
  else
    result_df = test_neutral_evolution( repetitions, numinputs, numinstructions_rng, numregisters_rng, max_steps, cartesian=cartesian  )
  end
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", default_funcs(numinputs))
      println(f,"# cartesian: ",cartesian)
      println(f,"# repetitions: ",repetitions)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, result_df, append=true, writeheader=true )
    end
  end
  result_df
end

function test_neutral_evolution( repetitions::Int64, p::Parameters, max_steps::Int64; 
      cartesian::Bool=:false, print_steps::Bool=:false )
  test_neutral_evolution( repetitions, p.numinputs, p.numinteriors, p.numlevelsback, max_steps, 
      cartesian=cartesian, print_steps=print_steps )
end

function test_neutral_evolution( repetitions::Int64, numinputs::Int64, numinstructions_rng::IntRange, numregisters_rng::IntRange,
      max_steps::Int64; cartesian::Bool=:false, print_steps::Bool=:false )
  df = DataFrame()
  df.numinputs = Int64[]
  if cartesian
    df.numgates = Int64[]
    df.levelsback = Int64[]
  else
    df.numinstructions = Int64[]
    df.numregisterss = Int64[]
  end
  df.steps = Float64[]
  df.fail_fract = Float64[]
  funcs = lin_funcs(numinputs)
  for ni in numinstructions_rng
    for nr in numregisters_rng
      println("nr: ",nr)
      p = Parameters(numinputs,1,ni,nr)
      sum_steps = 0
      sum_fail = 0
      for i = 1:repetitions
        gl = randgoal(p)
        if cartesian
          c = random_chromosome( p )
        else
          c = rand_lcircuit(p,funcs)
        end
        (c,step) = neutral_evolution( c, gl, max_steps, print_steps=print_steps )
        if repetitions > 1000 && i % 100 == 0
          println("i: ",i,"  gl: ",@sprintf("0x%x",gl[1]),"  ni: ",ni,"  nr: ",nr,"  step: ",step)
        elseif repetitions <= 1000
          println("i: ",i,"  gl: ",@sprintf("0x%x",gl[1]),"  ni: ",ni,"  nr: ",nr,"  step: ",step)
        end
        sum_steps += step
        sum_fail += step==max_steps ? 1 : 0
      end
      df_row = ( numinputs, ni, nr, sum_steps/repetitions, sum_fail/repetitions )
      push!(df,df_row)
    end
  end
  df
end

function average_compatible_dataframes( df1::DataFrame, df2::DataFrame )
  @assert names(df1) == names(df2)
  @assert size(df1) == size(df2)
  result_df = DataFrame()
  for i = 1:size(df2)[2]
    if df1[:,i] == df2[:,i]
      insertcols!(result_df,i,names(df1)[i]=>df1[:,i])
    else
      insertcols!(result_df,i,names(df1)[i]=>(df1[:,i]+df2[:,i])./2.0)
    end
  end
  result_df
end
