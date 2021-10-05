# Construct the phenotype network for a given parameter setting.
using HDF5, JLD
# Constructs the phenotype network matrix and outlist which is the counts
#   of genotypes that map to the corresponding phenotype before mutation
function construct_pheno_net_parallel( p::Parameters, nreps::Int64, numcircuits::Int64; 
    use_lincircuit::Bool=false, csvfile::String="", pheno_file::String="") 
  num_processes_per_processor = 4
  println("nreps: ",nreps,"  numcircuits: ",numcircuits)
  funcs = default_funcs(p.numinputs)
  if pheno_file == "" || !isfile(pheno_file)
    ph_net = create_empty_pheno_net( p )
    outlist = zeros(Int64,2^2^p.numinputs)
  elseif isfile(pheno_file)
    (outlist,ph_net)=jldopen(pheno_file,"r") do file
      (read( file, "outlist" ), read( file, "ph_net" ))
     end
  end
  #result_list = Tuple{Matrix{Int64},Vector{Int64}}[]
  pmap_list = nprocs() <= 1 ? [1] : collect(1:(num_processes_per_processor*(nprocs()-1)))
  denom = length(pmap_list)
  println("denom: ",denom,"  pmap_list: ",pmap_list)
  result_list = pmap(x->construct_pheno_net( p, div(nreps,denom), div(numcircuits,denom), use_lincircuit=use_lincircuit ), pmap_list )
  #println("result_list: ",result_list)
  ph_net = create_empty_pheno_net( p )
  outlist = zeros(Int64,2^2^p.numinputs)
  for (phn, outlst) in result_list
    ph_net .+= phn
    outlist .+= outlst
  end
  if pheno_file != ""
    save( pheno_file, "outlist", outlist, "ph_net", ph_net )
  end
  if csvfile != ""
    @assert csvfile[end-3:end] == ".csv"
    mdf = marginals_dataframe( ph_net, outlist, p )
    pn_df = pheno_net_df( ph_net, p )
    ph_file = csvfile[1:end-4] * "_phnet" * ".csv"
    write_csv_file( ph_file, pn_df, funcs, num_processes_per_processor )
    rc_file = csvfile[1:end-4] * "_rowcol" * ".csv"
    write_csv_file( rc_file, mdf, funcs, num_processes_per_processor )
  end
  (ph_net, outlist)
end

# Constructs the phenotype network matrix and outlist which is the counts
#   of genotypes that map to the corresponding phenotype before mutation
# Note:  pmap paralleism failed badly on 9/28/21
function construct_pheno_net( p::Parameters, nreps::Int64, numcircuits::Int64; 
    use_lincircuit::Bool=false) 
  println("nreps: ",nreps,"  numcircuits: ",numcircuits)
  funcs = default_funcs(p.numinputs)
  ph_net = create_empty_pheno_net( p )
  (outlist, ch_ints_list) = count_outputs( nreps, p, numcircuits, use_lincircuit=use_lincircuit )
  #(new_outlist, ch_ints_list) = count_outputs_parallel( nreps, p, numcircuits, use_lincircuit=use_lincircuit )
  #outlist .+= new_outlist
  #println("length(outlist): ",length(outlist))
  println("length(ch_ints_list): ",length(ch_ints_list))
  #println("ch_ints_list: ",ch_ints_list)
  #new_ph_net_list = pmap(ch_ints->process_ch_ints( ch_ints, p, funcs, use_lincircuit=use_lincircuit ), ch_ints_list ) # pmap commented out 10/4/21
  new_ph_net_list = map(ch_ints->process_ch_ints( ch_ints, p, funcs, use_lincircuit=use_lincircuit ), ch_ints_list )
  for new_ph_net in new_ph_net_list
    ph_net .+= new_ph_net
  end
  (ph_net, outlist)
end

# Helper function for construct_pheno_net()
# Note:  pmap paralleism failed badly on 9/28/21
function process_ch_ints( ch_ints_list::Vector{Int128} , p::Parameters, funcs::Vector{Func}; 
      use_lincircuit::Bool=false )
  #println("process_ch_ints: ")
  ph_net = create_empty_pheno_net( p )
  for ch_int in ch_ints_list
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
  ph_net
end  

function create_empty_pheno_net( p::Parameters; pheno_file::String="" )
  pheno_net = zeros( Int64, p.numoutputs*2^2^p.numinputs, p.numoutputs*2^2^p.numinputs )
  if pheno_file != ""
    save( pheno_file, "pheno_net", pheno_net )
  end
  pheno_net
end

function write_csv_file( csvfile::String, df::DataFrame, funcs::Vector{Func}, num_processes_per_processor::Int64 )
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ",funcs)
    println(f,"# nreps: ",nreps)
    println(f,"# numcircuits: ",numcircuits)
    println(f,"# num_processes_per_processor: ",num_processes_per_processor)
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

# Converts the ph_net matrix to a dataframe.  
# The first column is the goals
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

# returns the tuple (ph_net,outlist)
function read_jld_file( jld_file::String )
  (outlist,ph_net)=jldopen(jld_file,"r") do file
    (read( file, "ph_net" ), read( file, "outlist" )) 
  end
end  

# Approximates the Markov chain stationary distribution by post-multiplying init_distribution by the
#   transition matrix niters times.
# The Markov Chain transition matrix is computed by normalizing each row of the phenotype net matrix to sum to 1.
function markov_chain_stationary( niters::Int64, ph_net::Array{Int64,2}, init_distribution::AbstractArray )
  dim = size(ph_net)[1]
  println("dim: ",dim)
  @assert dim == size(ph_net)[2]
  if size(init_distribution)[1] == dim
    init_distribution = transpose(init_distribution)
  end
  @assert size(init_distribution)[1] == 1
  @assert size(init_distribution)[2] == dim
  nonzero_rows = transpose([ (sum(ph_net[i,:]) != 0 ? 1 : 0) for i = 1:dim ])
  init_distribution .*= nonzero_rows   # zero elements of init_distribution that correspond to zero rows of ph_net
  if sum(init_distribution) != 1.0
    init_distribution = init_distribution/sum(init_distribution)
  end
  println("init: ",init_distribution)
  T = zeros(Float64, dim, dim)  # Transition matrix
  for i = 1:dim
    T[i,:] =  sum(ph_net[i,:])>0 ? ph_net[i,:]/sum(ph_net[i,:]) : zeros(Float64,dim)
  end
  #println("T: ",T[1:20,1:20])
  for t = 1:niters
    init_distribution = init_distribution*T
    #println("init: ",init_distribution)
  end
  init_distribution
end

