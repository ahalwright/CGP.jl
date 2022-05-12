# Creates an array of the outputs created by nreps random chromosomes.
# Prints (shows) this array, showing elemets greater than a show_increment parameter, or less than or equal to  show_increment.
# Or shows the count of outputs satisfying these conditions.
# Since an array is used to store the outputs, numinputs for chromosomes is limited to 4.
# In order to handle larger numinputs, a hash table would need to be used---not implemented at this time.
# Example:
# julia> ol = count_outputs_parallel( 100, 2, 1, 4, 4 );
# julia> show_output_list( ol, 2 )  # shows outputs created more than 2 times
# output: 0X0:  24
# output: 0X1:  8
# output: 0X2:  3
# output: 0X3:  4
# output: 0X5:  11
# output: 0X6:  3
# output: 0X7:  6
# output: 0X8:  5
# output: 0X9:  3
# output: 0Xa:  13
# output: 0Xb:  7
# output: 0Xc:  6
# output: 0Xd:  4
# output: 0Xe:  8
# output: 0Xf:  18
# count of shown outputs: 15
# julia> show_outputs_list( ol, 2, show_small=true )  # shows outputs created at most twice
# output: 0X4:  2
# count of shown outputs: 1
# See notes/diargy9_9_20.txt for more examples.
using DataFrames
using Distributed
using Printf
using Dates
using CSV
export count_outputs, count_outputs_parallel, write_to_file, read_file, show_outputs_list, read_counts_files
export add_counts_to_dataframe, write_to_dataframe_file, circuit_complexities, run_circuit_complexities
export create_lincircuits_list, create_chromosome_ints_list, increment_circuit_ints_list! 

MyFunc = Main.CGP.MyFunc
# Return an output list of the number of times that an output was produced by randomly generating chromosomes with these parameters
#  This list is indexed over phenotypes.
# If output_complex==true the returned outlist is a list of pairs, where the first element of the pair is the number of times the
#   output is produced, and the second is the average of the Tononi complexities of the genotypes that produce that output.
# Not used by construct_pheno_net() in PhenotypeNetwork.jl since nesting calls to pmap() doesn't work (stack overflow error)
function count_outputs_parallel( nreps::Int64, p::Parameters, numcircuits::Int64, funcs::Vector{Func}=Func[]; 
    csvfile::String="", use_lincircuit::Bool=:false, output_complex::Bool=false )
  count_outputs_parallel( nreps::Int64, p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, numcircuits, funcs,
    csvfile=csvfile, use_lincircuit=use_lincircuit, output_complex=output_complex )
end

# Not used by construct_pheno_net() in PhenotypeNetwork.jl since nesting calls to pmap() doesn't work (stack overflow error)
function count_outputs_parallel( nreps::Int64, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64, numcircuits::Int64, 
    funcs::Vector{Func}=Func[]; csvfile::String="", use_lincircuit::Bool=:false, output_complex::Bool=false ) 
  #=
  if use_lincircuit && output_complex
    println("Warning: output_complex reset to false when use_lincircuit==true")
    output_complex=false
  end
  =#
  p = Parameters( numinputs, numoutputs, numinteriors, numlevelsback )
  if length(funcs) == 0
    funcs = use_lincircuit ? lin_funcs(numinputs) : default_funcs(numinputs)
  end
  print_parameters( p )
  println("nprocs: ",nprocs())
  #n_procs=nprocs()
  if nprocs() > 1
    nreps_p = Int(round(nreps/(nprocs()-1)))
  else
    nreps_p = nreps 
  end
  println("nreps_p: ",@sprintf("%.1e",nreps_p))
  println("numcircuits: ",numcircuits)
  #println("csvfile: ",csvfile) 
  count_out_funct(x) = count_outputs( nreps_p, numinputs, numoutputs, numinteriors, numlevelsback, numcircuits, 
      use_lincircuit=use_lincircuit, output_complex=output_complex )
  result =  pmap( x->count_out_funct(x), collect(1:nprocs()))
  #result =  map( x->count_out_funct(x), collect(1:nprocs()))
  println("len(result): ",length(result))
  #println("result[1][1]: ",result[1][1])
  #println("result[1][2]: ",result[1][2])
  outlist = reduce(+,map(x->x[1],result))   # Works for output_complex and !output_complex
  if output_complex   # Average complexities
    for i in 1:length(outlist)
      pair = outlist[i]
      new_pair = (pair[1], pair[2]/pair[1] )
      outlist[i] = new_pair
    end
  end
  circ_ints_list = result[1][2]
  for i = 2:length(result)
    vcat_arrays!(circ_ints_list,result[i][2])
  end
  if length(csvfile) > 0
    write_to_dataframe_file( p, outlist, circ_ints_list, funcs, numcircuits, nreps, csvfile=csvfile )
  end
  #println("outlist: ",outlist)
  #println("circ_ints_list: ",circ_ints_list)
  (outlist,circ_ints_list)
end

# Return an output list of the number of times that an output (phenotype)  was produced by randomly generating genotypes with these parameters
#  This list is indexed over phenotypes.
# If output_complex==true the returned outlist is a list of pairs, where the first element of the pair is the number of times the
#   output is produced, and the second is the sum of the Tononi complexities of the genotypes that produce that output.
function count_outputs( nreps::Int64, p::Parameters, numcircuits::Int64=0; use_lincircuit::Bool=:false, output_complex::Bool=false )
  count_outputs( nreps::Int64, p.numinputs, p.numoutputs, p.numinteriors, p.numlevelsback, numcircuits, 
    use_lincircuit=use_lincircuit, output_complex=output_complex )
end

# Return an output list of the number of times that an output (phenotype) was produced by randomly generating genotypes with these parameters
function count_outputs( nreps::Int64, numinputs::Int64, numoutputs::Int64, numinteriors::Int64, numlevelsback::Int64, numcircuits::Int64=0;
    use_lincircuit::Bool=:false, output_complex::Bool=false )
  p = Parameters( numinputs=numinputs, numoutputs=numoutputs, numinteriors=numinteriors, numlevelsback=numlevelsback ) 
  funcs = use_lincircuit ? lin_funcs(numinputs) : default_funcs(numinputs)
  if output_complex
    outlist = fill( (0,0.0), numoutputs*2^2^numinputs )
  else
    #outlist = fill( Int64(0), numoutputs*2^2^numinputs )
    outlist = fill( UInt128(0), numoutputs*2^2^numinputs )
  end
   circuit_ints_list = [ Int128[] for _ = 1:p.numoutputs*2^2^p.numinputs ]
  for _ = 1:nreps
    if use_lincircuit
      c = rand_lcircuit( p, funcs )
      c_int = instruction_ints_to_circuit_int( instruction_vects_to_instruction_ints( c, funcs ), p, funcs )
      @assert c_int >= 0
      #=
      if c_int < 0 
        println(" c_int: ",c_int,"  c: ",c)
      end
      =#
    else
      c = random_chromosome( p, funcs )
      c_int = chromosome_to_int( c, funcs)
    end
    output = output_values( c )
    increment_circuit_ints_list!( circuit_ints_list, output, c_int, numcircuits, p, funcs ) 
    if output_complex
      #println("outlist: ",outlist)
      complexity = use_lincircuit ? lincomplexity(c,funcs) : complexity5(c)
      #println("outlist[concatenate_outputs(output,numinputs)+1][1]: ",outlist[concatenate_outputs(output,numinputs)+1][1])
      pair = outlist[concatenate_outputs(output,numinputs)+1]
      new_pair = (pair[1]+1,pair[2]+complexity)
      outlist[concatenate_outputs(output,numinputs)+1] = new_pair
    else
      outlist[concatenate_outputs(output,numinputs)+1] += 1
    end
  end
  (outlist,circuit_ints_list)
end

import Base.:+
# Define + on outlist tuples
function +(t1::Tuple{Int64,Float64},t2::Tuple{Int64,Float64})
  (t1[1]+t2[1],t1[2]+t2[2])
end 

# increments the circuit_list corresponding to output if it contains less than numcircuits circuits
function increment_circuit_ints_list!( circuit_ints_list::Vector{Vector{Int128}}, output::Goal, c_int::Int128, 
      numcircuits::Int64, p::Parameters, funcs::Vector{Func} )
  index = concatenate_outputs(output,p.numinputs)+1
  if length(circuit_ints_list[index]) >= numcircuits
    return
  end
  # c_ints = map(lc->instruction_vect_to_instruction_int(lc,numregisters,numinputs,funcs), circ.circuit_vects)
  #println("index: ",index,"  output: ",output,"  execute: ", execute_lcircuit( c_ints, numregisters, numinputs, funcs ))
  push!(circuit_ints_list[index], c_int )
end

# increments the chrome_ints_list corresponding to output if it contains less than numcircuits chromosome ints
# Not used by function count_outputs()
function increment_chrome_ints_list!( chrome_ints_list::Vector{Vector{Int128}}, output::Goal, ch::Chromosome, numcircuits::Int64,
    p::Parameters, funcs::Vector{Func} )
  index = concatenate_outputs(output,p.numinputs)+1
  if length(chrome_ints_list[index]) >= numcircuits
    return
  end
  ch_int = chromosome_to_int(ch)
  push!(chrome_ints_list[index], ch_int )
end

# If numoutputs > 1, combines the outputs into a single unsigned integer of type MyFunc
# If numoutputs == 1, extract the one component
function concatenate_outputs( output::Goal, numinputs::Int64 )
  if length(output) == 1
    return MyFunc(output[1])
  else
    result = MyFunc(0)
    for i = 1:length(output)
      result <<= 2^numinputs
      result |= output[i]
    end
    result
  end
end

# Example:  mycat([[[3,2],[5,9]],[[8,4],[6,7]]])
#    [[3,2],[5,9],[8,4],[6,7]]
function mycat( lsts::Vector{Vector{Vector{Int64}}} ) 
  result = Vector{Int64}[]
  for lst in lsts
    if length(lst) > 0
      result = vcat(result,lst)
    end
  end
  result
end
function mycat( lsts::Vector{Vector{Int128}} ) 
  result = Vector{Int64}[]
  for lst in lsts
    if length(lst) > 0
      result = vcat(result,lst)
    end
  end
  result
end

# Concatenates (using vcat()) the components of lst1 and lst2 with the result replacing lst1
# Example:  
#  vcat_arrays!([Int128[5,9,4],Int128[3,1,7]],[Int128[8,6],Int128[2]]) 
#    [Int128[5,9,4,8,6],Int128[3,1,7,2]]
function vcat_arrays!( lst1::Vector{Vector{Int128}}, lst2::Vector{Vector{Int128}} )
  @assert length(lst1) == length(lst2)
  for i = 1:length(lst1)
    lst1[i] = vcat(lst1[i],lst2[i])
  end
  lst1
end
#=
function write_to_dataframe_file( p::Parameters, outputs_list::Vector{UInt128}, funcs::Vector{Func}, numcircuits::Int64=0, nreps::Int64=0; csvfile::String="" )
  df = DataFrame()
  df.:goals = [ @sprintf("0x%04x",g) for g = 0:(2^2^p.numinputs-1) ]
  #println("len goals: ",length(df.:goals))
  sym = Symbol("ints","$(p.numinteriors)","_","$(p.numlevelsback)") 
  df[!,sym] = outputs_list
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end) 
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )     
    print_parameters(f,p,comment=true)
    if numcircuits > 0
      println(f,"# numcircuits: ",numcircuits)
    end
    if nreps > 0
      spr = @sprintf("%.2e",nreps)
      println(f,"# nreps: $(spr)")
    end
    println(f,"# funcs: ",funcs)
    CSV.write( f, df, append=true, writeheader=true )
  end
end

#function write_to_dataframe_file( p::Parameters, outputs_list::Vector{Int64}, funcs::Vector{Func}; csvfile::String="" )
function write_to_dataframe_file( p::Parameters, outputs_list::Vector{UInt128}, funcs::Vector{Func}, numcircuits::Int64=0, nreps::Int64=0; csvfile::String="" )
  df = DataFrame()
  df.:goals = [ @sprintf("0x%04x",g) for g = 0:(2^2^p.numinputs-1) ]
  df.:igoals = [ @sprintf("%d",g) for g = 0:(2^2^p.numinputs-1) ]
  #println("len goals: ",length(df.:goals))
  sym = Symbol("ints","$(p.numinteriors)","_","$(p.numlevelsback)") 
  df[!,sym] = outputs_list
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end) 
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )     
    print_parameters(f,p,comment=true)
    if numcircuits > 0
      println(f,"# numcircuits: ",numcircuits)
    end
    if nreps > 0
      spr = @sprintf("%.2e",nreps)
      println(f,"# nreps: $(spr)")
    end
    println(f,"# funcs: ",funcs)
    CSV.write( f, df, append=true, writeheader=true )
  end
end
=#

function write_to_dataframe_file( p::Parameters, outputs_list::Vector{UInt128}, circuits_list::Vector{Vector{Vector{Int64}}}, funcs::Vector{Func}, numcircuits::Int64=0, nreps::Int64=0;
     csvfile::String="" )
  df = DataFrame()
  df.:goals = [ @sprintf("0x%04x",g) for g = 0:(2^2^p.numinputs-1) ]
  #println("len goals: ",length(df.:goals))
  sym = Symbol("ints","$(p.numinteriors)","_","$(p.numlevelsback)") 
  df[!,sym] = outputs_list
  df.:circuits_list = circuits_list
  if length(csvfile) > 0
    write_df_to_csv( df, p, funcs, csvfile, numcircuits=numcircuits, nreps=nreps )
  end
  df 
end

function write_to_dataframe_file( p::Parameters, outputs_list::Vector{UInt128}, circuits_list::Vector{Vector{Int128}}, funcs::Vector{Func}, numcircuits::Int64=0, nreps::Int64=0; 
    csvfile::String="" )
  df = DataFrame()
  df.:goals = [ @sprintf("0x%04x",g) for g = 0:(2^2^p.numinputs-1) ]
  #println("len goals: ",length(df.:goals))
  sym = Symbol("ints","$(p.numinteriors)","_","$(p.numlevelsback)") 
  df[!,sym] = outputs_list
  df.:circuits_list = circuits_list
  if length(csvfile) > 0
    write_df_to_csv( df, p, funcs, csvfile, numcircuits=numcircuits, nreps=nreps )
  end
  df 
end

function write_to_dataframe_file( p::Parameters, outputs_list::Vector{Tuple{Int64,Float64}}, funcs::Vector{Func}, numcircuits::Int64=0, nreps::Int64=0; csvfile::String="" )
  df = DataFrame()
  df.:goals = [ @sprintf("0x%04x",g) for g = 0:(2^2^p.numinputs-1) ]
  sym = Symbol("ints","$(p.numinteriors)","_","$(p.numlevelsback)") 
  df[!,sym] = map(x->x[1],outputs_list)
  df.:complexity = map( x->x[2], outputs_list )
  #df.:circuits_list = circuits_list
  if length(csvfile) > 0
    write_df_to_csv( df, p, funcs, csvfile, numcircuits=numcircuits, nreps=nreps )
  end
  df
end

#= Moved to Utilities.jl on 4/20/22
function write_df_to_csv( df::DataFrame, p::Parameters, funcs::Vector{Func}, csvfile::String )
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end) 
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )     
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ",funcs)
    CSV.write( f, df, append=true, writeheader=true )
  end
end
=#

# print the elements of the output list that are greater than increment
# if show_small==true, print the elements of the output list that are less than or equal to  increment
function show_outputs_list( outputs_list::Vector{MyFunc}, show_increment::Int64=0; show_small::Bool=false, count_only::Bool=false )
  count = 0
  for i = 0:length(outputs_list)-1
    if !show_small
      if outputs_list[i+1] > show_increment 
        if !count_only
          @printf("output: 0X04%x:  %d\n",i,outputs_list[i+1])
        end
        count += 1
      end
    else
      if outputs_list[i+1] <= show_increment 
        if !count_only
          @printf("output: 0X04%x:  %d\n",i,outputs_list[i+1])
        end
        count += 1
      end
    end
  end
  println("count of shown elements: ",count)
end

# Write output counts list to a file.  
# Optional list of comments which are written to a file preceded by a # character.
# Example:  write_to_file( ol, "../data/9_8/test.csv", "10^6 reps Raman funcs", "3 inputs, 1 output, 7 interiors, 4 levsback",hex=true )
function write_to_file( outputs_list::Vector{MyFunc}, filename::String, comments::String... ; 
      hex::Bool=false )
  open( filename, "w" ) do f
    for s in comments
      println(f,"# ",s)
    end
    for i = 1:length(outputs_list)
      if hex
        @printf(f,"0x%04x\n",outputs_list[i])
      else
        println(f, outputs_list[i] )
      end
    end
  end
end

# Read a counts file written by write_to_file(). 
# Returns a list whose elements are of type out_type.
# Options for out_type include Int64, UInt64, UInt32,
# Lines starting with "#" are comments
# Didn't work 10/11/20
function read_file( filename::String, out_type::Type )
  result = out_type[]
  i = 0
  for row in CSV.File(filename, header=false, comment="#" )
    #println("i: ",i,"  ",row.Column1)
    push!(result,parse(out_type,row.Column1))
    i += 1
  end
  result
end

# Read multiple counts files where each counts file has one base 10 integer per line.
# Returns a 2-tuple where the first element is a list of goals written in hex form,
#    and the second is a vector of counts lists, one for each file read.
function read_counts_files(filename::String...)
  my_goals = Vector{String}[]
  my_counts = Vector{Int64}[]
  i=1
  for fname in filename
    open(fname) do f
      lines = readlines(f)
      if i > 1
        @assert length(lines) == length(my_goals[i-1])
      end
      filter!(x->x[1]!='#',lines)
      goals = [ @sprintf("0x%04x",i) for i = 0:(length(lines)-1) ]
      if i > 1
        @assert goals == my_goals[i-1]
      end
      push!(my_goals,goals)
      counts = [ parse(Int64, lines[i+1] ) for i = 0:(length(lines)-1)]
      push!(my_counts,counts)
    end
  end
  (my_goals[1],my_counts)
end

# Add multiple vectors of numbers as additional columns to dataframe df.
# The names of the columns are given in the names vector.
# The filename to read the counts vector is filename.
# The column file has one entry which must be a base 10 integer.
# The dataframe can be read using function read_dataframe() in Analyze.jl
function add_counts_to_dataframe( df::DataFrame, names::Vector{Symbol}, filename::String... )
  (goals,counts) = read_counts_files( filename... )
  df.goal = goals
  i = 1
  for nm in names
    insertcols!(df, size(df)[2]+1, nm=>counts[i] )
    i += 1
  end
  df
end

# The counts dataframe is read from the file counts_filename
# Assumes that the counts dataframe has a field "goals" of type String.  Examples of this field:  "0xe","0x10","0xffee".
# Assumes that the df dataframe has a filed "goal" of type String.  Example of this field:  "UInt16[0x09b5]".
function add_counts_to_dataframe( df::DataFrame, counts_filename::String, counts_field::Symbol )
  cdf = read_dataframe(counts_filename)
  # eval(Meta.parse(df.goal[i]))  converts the string df.goal[i] to Julia Vector.  Example:  UInt16[0x09b5]
  # eval(Meta.parse(df.goal[i]))[1]  extracts the sole element of this vector
  # @sprintf("0x%04x",eval(Meta.parse(df.goal[i]))[1])   converts to a string of the same format as in the goals field of cdf.
  counts = [cdf[cdf.goal.==@sprintf("0x%04x",eval(Meta.parse(df.goal[i]))[1]),counts_field][1] for i = 1:size(df)[1]]
  #println("counts: ",counts)
  insertcols!(df, size(df)[2]+1, counts_field=>counts )
  df
end 

function merge_count_dataframes( df1::DataFrame, df2::DataFrame )
  df = DataFrame()
  if  "goals" in names(df1) && "goals" in names(df2)
    df.goals = df1.goals
  end
  if "ints6_2" in names(df1) && "ints6_2" in names(df2)
    df.ints6_2 = df1.ints6_2 + df2.ints6_2
  end
  if "circuits_list" in names(df1) && "circuits_list" in names(df2)
    df1.cl = map(x->eval(Meta.parse(x)),df1.circuits_list)
    df2.cl = map(x->eval(Meta.parse(x)),df2.circuits_list)
    df.circuits_list = map((x,y)->vcat(x,y),df1.cl,df2.cl)
  end
  df
end

# Create a pair (goal_list,complexity_list) where a pair (goal_list[i], complexity_list[i]) corresponds to a random circuit.
function circuit_complexities( p::Parameters, num_circuits::Int64 )
  println("circuit_complexities: num_circuits: ",num_circuits)
  funcs = default_funcs(p.numinputs)
  complexity_list = zeros(Float64,num_circuits)
  goal_list = Goal[]
  for i = 1:num_circuits
    c = random_chromosome( p, funcs )
    push!( goal_list, output_values(c) )
    complexity = complexity5(c)
    complexity_list[i] = complexity >= 0.0 ? complexity : 0.0
  end
  ( goal_list, complexity_list )
  #=
  df = DataFrame()
  df.goal = goal_list
  df.complexity = complexity_list
  sort!( df, [:complexity], rev=true )
  df
  =#
end

function run_circuit_complexities( p::Parameters, num_circuits::Int64; csvfile::String="" )
  num_circuits_per_proc = Int(trunc(num_circuits/(nprocs()-1)))
  goal_complexity_pairs = pmap(x->circuit_complexities( p, num_circuits_per_proc), collect(1:(nprocs()-1)) )
  goals = Goal[]
  complexities = Float64[]
  for gc in goal_complexity_pairs
    goals = vcat( goals, gc[1] )
    complexities = vcat( complexities, gc[2] )
  end
  df = DataFrame()
  df.goals = goals
  df.complexities = complexities
  if length(csvfile) > 0
    open( csvfile, "w" ) do f 
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", Main.CGP.default_funcs(p.numinputs))
      print_parameters(f,p,comment=true)
      println(f,"# num_circuits: ",num_circuits) 
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

#=
function check_geno_exist_for_phenotypes( p::Parameters, gl::GoalList, max_evolve_tries::Int64, max_ev_steps::Int64, outfile::String )
  funcs = default_funcs(p.numinputs)
  found_list = Vector{UInt16}[]
  not_found_list = Vector{UInt16}[]
  for g in gl
    i = 1
    while i < max_evolve_tries
      println("i: ",i)
      c = random_chromosome(p,funcs)
      (c,step,worse,same,better,output,matched_goals,matched_goals_list) =
          mut_evolve( c, [g], funcs, max_ev_steps, print_steps=false )
      if step < max_ev_steps
        break
      end
      i += 1
    end
    if i < max_evolve_tries
      @printf("circuit found for goal: [0x%04x]\n",g[1])
      push!(found_list,g)
    else
      @printf("circuit not found for goal: [0x%04x]\n",g[1])
      push!(not_found_list,g)
    end
  end  
  open( outfile, "w" ) do f
    print_parameters(f,p)
    println(f,"max_evolve_tries: ",max_evolve_tries)
    println(f,"max_ev_steps: ",max_ev_steps)
    println(f,"found_list")
    println(f,found_list)
    println(f,"not found_list")
    println(f,not_found_list)
  end
end
=#

