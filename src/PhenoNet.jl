# Construct the phenotype network for a given parameter setting.
using HDF5, JLD
using Base.Threads
# The phenotype network has phenotypes as nodes and mutations between phenotypes as edges.
# Edges are weighted by the number of mutations from the source phenotype to the destination phenotype.
# However, mutations for an edge from a phenotype are limited in number by the parameter numcircuits.
# This is because we are not so interested in common phenotypes, so if we under-weight edges from common phenotypes,
#   this won't skew are results.
# The outlist is the number of copies of each phenotype that are sampled before mutation.  It is not limited by the numcircuits parameter.
# Constructs the phenotype network matrix and outlist.
function construct_pheno_net_parallel( p::Parameters, nreps::Int64, numcircuits::Int64; 
    use_lincircuit::Bool=false, csvfile::String="", pheno_file::String="") 
  global sync_flag_outlist = 0
  global sync_flag_matrix = 0
  nthreads = Threads.nthreads()
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
  outlist = zeros(Int64,2^2^p.numinputs)
  cndition_outlist = Threads.Condition()
  cndition_matrix = Threads.Condition()
  Threads.@threads for i = 1:nthreads
    increment_pheno_net!( ph_net, outlist, p, cndition_outlist, cndition_matrix, 
        div(nreps,nthreads), div(numcircuits,nthreads), use_lincircuit=use_lincircuit )
  end
  if pheno_file != ""
    save( pheno_file, "outlist", outlist, "ph_net", ph_net )
  end
  if csvfile != ""
    @assert csvfile[end-3:end] == ".csv"
    mdf = marginals_dataframe( ph_net, outlist, p )
    #println("mdf: ",mdf)
    pn_df = pheno_net_df( ph_net, p )
    ph_file = csvfile[1:end-4] * "_phnet" * ".csv"
    write_csv_file( ph_file, pn_df, funcs, nthreads )
    rc_file = csvfile[1:end-4] * "_rowcol" * ".csv"
    write_csv_file( rc_file, mdf, funcs, nthreads )
  end
  (ph_net, outlist)
end

# Increments the phenotype network matrix and outlist which is the counts
#   of genotypes that map to the corresponding phenotype before mutation
# Modifies both ph_net and outlist while their conditions are nlocked 
function increment_pheno_net!( ph_net::Matrix{Int64}, outlist::Vector{Int64}, p::Parameters, 
      cndition_matrix::Base.GenericCondition{ReentrantLock}, cndition_outlist::Base.GenericCondition{ReentrantLock}, 
      nreps::Int64, numcircuits::Int64; use_lincircuit::Bool=false) 
  global sync_flag_outlist
  global sync_flag_matrix
  println("nreps: ",nreps,"  numcircuits: ",numcircuits)
  funcs = default_funcs(p.numinputs)
  # The count_outputs() function is in RecordOutputs.jl
  (new_outlist, ch_ints_list) = count_outputs( nreps, p, numcircuits, use_lincircuit=use_lincircuit )
  lock(cndition_outlist)
  try
    while sync_flag_outlist != 0
      wait(cndition_outlist)
    end
  finally
    sync_flag_outlist = sync_flag_outlist + 1
    outlist .+= new_outlist
    sync_flag_outlist = sync_flag_outlist - 1
    unlock(cndition_outlist)
  end  
  println("length(ch_ints_list): ",length(ch_ints_list))
  for ch_ints in ch_ints_list
    process_ch_ints!( ph_net, cndition_matrix, ch_ints, p, funcs, use_lincircuit=use_lincircuit )
  end
  nothing
end

# Helper function for construct_pheno_net()
# Modifies ph_net while cndition_matrix is unlocked.
function process_ch_ints!( ph_net::Matrix{Int64}, cndition_matrix::Base.GenericCondition{ReentrantLock}, ch_ints_list::Vector{Int128} ,
      p::Parameters, funcs::Vector{Func}; use_lincircuit::Bool=false )
  global sync_flag_matrix
  for ch_int in ch_ints_list
    if use_lincircuit
      ch = circuit_int_to_circuit( ch_int, p, funcs )
    else
      ch = int_to_chromosome( ch_int, p, funcs )
    end
    src = output_values(ch)[1]   # assumes 1 output    
    dests = map(x->x[1], mutate_all( ch, funcs ))
    lock(cndition_matrix)
    try
      while sync_flag_matrix != 0
        wait(cndition_matrix)
      end
    finally
      sync_flag_matrix = sync_flag_matrix + 1
      for dest in dests
        ph_net[src+1,dest+1] += 1
      end
      sync_flag_matrix = sync_flag_matrix - 1
      unlock(cndition_matrix)
    end
  end  
  nothing
end  

function create_empty_pheno_net( p::Parameters; pheno_file::String="" )
  pheno_net = zeros( Int64, p.numoutputs*2^2^p.numinputs, p.numoutputs*2^2^p.numinputs )
  if pheno_file != ""
    save( pheno_file, "pheno_net", pheno_net )
  end
  pheno_net
end

function write_csv_file( csvfile::String, df::DataFrame, funcs::Vector{Func}, nthreads::Int64 )
  open( csvfile, "w" ) do f
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    println(f,"# date and time: ",Dates.now())
    println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
    print_parameters(f,p,comment=true)
    println(f,"# funcs: ",funcs)
    println(f,"# use_lincircuit: ",use_lincircuit)
    println(f,"# nreps: ",nreps)
    println(f,"# numcircuits: ",numcircuits)
    println(f,"# nthreads: ",nthreads)
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
  result = evolvability_from_ph_net(ph_net,p)
  mdf.s_evolvability = result[1]
  mdf.d_evolvability = result[2]
  mdf.robustness = result[3]
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

function evolvability_from_ph_net( ph_net::Array{Int64,2}, p::Parameters )
  sigma( k::Int64 ) = ( k == 0 ? 0 : 1 )
  s_evolvability = [ sum( ph_net[i,:] ) - ph_net[i,i]  for i = 1:p.numoutputs*2^2^p.numinputs ]
  d_evolvability = [ sum( map(x->sigma(x), ph_net[i,:])) - sigma(ph_net[i,i])  for i = 1:p.numoutputs*2^2^p.numinputs ]
  robustness = [ ph_net[i,i]/sum( ph_net[i,:] )  for i = 1:p.numoutputs*2^2^p.numinputs ]
  (s_evolvability,d_evolvability,robustness)
end

function read_ph_net_from_csv( csvfile::String, p::Parameters )
  df = read_dataframe( csvfile )
  pheno_net = zeros( Int64, p.numoutputs*2^2^p.numinputs, p.numoutputs*2^2^p.numinputs )
  for i = 1:p.numoutputs*2^2^p.numinputs
    pheno_net[:,i] = df[:,i+1]
  end
  pheno_net
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

