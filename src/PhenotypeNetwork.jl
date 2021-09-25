# Construct the phenotype network for a given parameter setting.
using HDF5, JLD

function create_empty_pheno_net( p::Parameters; pheno_file::String="" )
  pheno_net = zeros( Int64, p.numoutputs*2^2^p.numinputs, p.numoutputs*2^2^p.numinputs )
  if pheno_file != ""
    save( pheno_file, "pheno_net", pheno_net )
  end
  pheno_net
end

function construct_pheno_net( p::Parameters, nreps::Int64, numcircuits::Int64; 
    use_lincircuit::Bool=false, csvfile::String="", pheno_file::String="") 
  funcs = default_funcs(p.numinputs)
  if pheno_file == "" || !isfile(pheno_file)
    ph_net = create_empty_pheno_net( p )
    outlist = zeros(Int64,2^2^p.numinputs)
  elseif isfile(pheno_file)
    #=
    ph_net = jldopen(pheno_file,"r") do file
      read( file, "pheno_net" )
    end 
    =#
    (outlist,ph_net)=jldopen(pheno_file,"r") do file
      (read( file, "outlist" ), read( file, "ph_net" ))
     end
  end
  #(outlist, ch_ints_list) = count_outputs( nreps, p, numcircuits, use_lincircuit=use_lincircuit )
  (new_outlist, ch_ints_list) = count_outputs_parallel( nreps, p, numcircuits, use_lincircuit=use_lincircuit )
  outlist .+= new_outlist
  println("length(outlist): ",length(outlist))
  println("length(ch_ints_list): ",length(ch_ints_list))
  #println("ch_ints_list: ",ch_ints_list)
  for ch_ints in ch_ints_list
    for ch_int in ch_ints
      #println("ch_int: ",ch_int)
      if use_lincircuit
        #ch_vects = instruction_ints_to_instruction_vects( circuit_int_to_instruction_ints( ch_int, p, funcs ), p, funcs )  
        #ch = LinCircuit( ch_vects, p ) 
        ch = circuit_int_to_circuit( ch_int, p, funcs )
      else
        ch = int_to_chromosome( ch_int, p, funcs )
        #print_circuit(ch)
      end
      src = output_values(ch)[1]   # assumes 1 output    
      #println("ch_int: ",ch_int,"  src: ",src)
      dests = map(x->x[1], mutate_all( ch, funcs ))
      #println("dests: ",dests)
      for dest in dests
        ph_net[src+1,dest+1] += 1
        #println("dest: ",dest,"  ph_net[src+1,dest+1]: ",ph_net[src+1,dest+1])
      end
    end
  end
  if pheno_file != ""
    save( pheno_file, "outlist", outlist, "ph_net", ph_net )
  end
  if csvfile != ""
    @assert csvfile[end-3:end] == ".csv"
    mdf = marginals_dataframe( ph_net, outlist, p )
    pn_df = pheno_net_df( ph_net, p )
    ph_file = csvfile[1:end-4] * "_phnet" * ".csv"
    write_csv_file( ph_file, pn_df, funcs )
    rc_file = csvfile[1:end-4] * "_rowcol" * ".csv"
    write_csv_file( rc_file, mdf, funcs )
  end
  (ph_net, outlist)
end

function write_csv_file( csvfile::String, df::DataFrame, funcs::Vector{Func} )
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ",funcs)
    println(f,"# nreps: ",nreps)
    println(f,"# numcircuits: ",numcircuits)
    CSV.write( f, df, append=true, writeheader=true )
  end
end

function marginals_dataframe( ph_net::Array{Int64,2}, outlist::Vector{Int64}, p::Parameters )
  @assert p.numoutputs == 1
  mdf = DataFrame()
  mdf.indices = [ i for i = 0:2^2^p.numinputs-1 ]  # Decimal indices
  mdf.goal = [ @sprintf("0x%04x",i) for i = 0:2^2^p.numinputs-1 ]  # MyInt indices
  mdf.outlist = outlist
  mdf.row_sums = [ sum( ph_net[i,:] ) for i = 1:2^2^p.numinputs ]
  mdf.col_sums = [ sum( ph_net[:,i] ) for i = 1:2^2^p.numinputs ] 
  mdf
end

function pheno_net_df( ph_net::Array{Int64,2}, p::Parameters ) 
  phdf = DataFrame()
  goals = [ @sprintf("0x%04x",i) for i = 0:2^2^p.numinputs-1 ]
  phdf.goal = goals
  for i = 1:2^2^p.numinputs
    colname = Symbol(goals[i])
    DataFrames.insertcols!(phdf, i+1, colname=>ph_net[:,i] )
  end
  phdf
end

function read_jld_file( jld_file::String )
  (outlist,ph_net)=jldopen(pheno_file,"r") do file
    (read( file, "ph_net" ), read( file, "outlist" )) 
  end
end  
