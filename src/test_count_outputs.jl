# Simple replication of count_outputs() in RecordOutputs.jl

function count_phenotypes_parallel( p::Parameters, funcs::Vector{Func}, nsamples::Int64; use_lincircuit::Bool=false, complexity::Bool=false, csvfile::String="" )
  numinputs = p.numinputs
  nsamples_iter = nprocs() > 1 ? Int64(ceil(nsamples/(nprocs()-1))) : nsamples
  println("nsamples_iter: ",nsamples_iter)
  iters_range = (nprocs()>1) ? (1:(nprocs()-1)) : (1:1)
  if complexity
    ph_cmplx_list = pmap( x->count_phenos( p, funcs, nsamples_iter, use_lincircuit=use_lincircuit, complexity=complexity ), iters_range )
    #pairs_list = map( x->count_phenos( p, funcs, nsamples_iter, use_lincircuit=use_lincircuit, complexity=complexity ), iters_range )
    phenos_list = map( x->x[1], ph_cmplx_list )
    complexities_list = map( x->x[2], ph_cmplx_list )
    phenos = reduce(+,phenos_list)
    cmblist = map(i->phenos_list[i].*complexities_list[i], 1:length(phenos_list))
    complexities = reduce(+,cmblist) ./ reduce(+,phenos_list)
    complexities = map(x-> isnan(x) ? 0 : x, complexities )
  else
    phenos_list = pmap( x->count_phenos( p, funcs, nsamples_iter, use_lincircuit=use_lincircuit ), iters_range )
    #phenos_list = map( x->count_phenos( p, funcs, nsamples_iter, use_lincircuit=use_lincircuit ), iters_range )
    phenos = reduce(+,phenos_list)
  end
  println("pmap reduce finished")
  df = DataFrame()
  goals = map(MyInt, 0:(2^2^numinputs-1) )
  df.goal = map(x->@sprintf("0x%04x",x),goals)
  df.counts = phenos
  if complexity
    df.complexity = complexities
  end
  println("size(df): ",size(df))
  if length(csvfile) > 0
    write_df_to_csv( df, p, funcs, csvfile, nsamples=nsamples )
  end
  df
end

function weighted_average_complexity( phcount1::Int64, cmplex1::Float64, phcount2::Int64, cmplex2::Float64 )
  (cmplex1*phcount1 + cmplex2*phcount2)/(phcount1+phcount2)
end

function accum_phenos_complexities( phenos_list::Vector{Vector{Int64}}, complexities_list::Vector{Vector{Float64}} )
  phenos = reduce(+,phlist)
  cmblist = map(i->phlist[i].*cplist[i], 1:length(phlist))
  complexities = reduce(+,cmblist) ./ reduce(+,phlist)
end

function count_phenotypes(  p::Parameters, funcs::Vector{Func}, nsamples::Int64; use_lincircuit::Bool=false, complexity::Bool=false, csvfile::String="" )
  numinputs = p.numinputs
  if complexity
    (phenos,complexities) = count_phenos( p, funcs, nsamples, use_lincircuit=use_lincircuit, complexity=complexity ) 
  else
    phenos = count_phenos( p, funcs, nsamples, use_lincircuit=use_lincircuit )
  end
  df = DataFrame()
  goals = map(MyInt, 0:(2^2^numinputs-1) )
  df.goal = map(x->@sprintf("0x%04x",x),goals)
  df.counts = phenos
  if complexity
    df.complexity = complexities
  end
  if length(csvfile) > 0
    write_df_to_csv( df, p, funcs, csvfile, ngens=nsamples )
  end
  df
end

function count_phenos(  p::Parameters, funcs::Vector{Func}, nsamples::Int64; use_lincircuit::Bool=false, complexity::Bool=false )
  #Random.seed!(1)   # use to test the equivalence of count_phenotypes() and count_phenotypes_parallel()
  numinputs = p.numinputs
  phenos = zeros(Int64,2^2^numinputs)
  avg_complexity = zeros(Float64,2^2^numinputs)
  for i = 1:nsamples
    if use_lincircuit
      c = rand_lcircuit( p, funcs )
    else
      c = random_chromosome( p, funcs )
    end
    ph = output_values( c )[1]
    phenos[ph+1] += 1
    if complexity
      cmplx = use_lincircuit ? lincomplexity(c,funcs) : complexity5(c)
      ph_count = phenos[ph+1]
      avg_complexity[ph+1] = ((ph_count-1)*avg_complexity[ph+1] + cmplx)/ph_count
    end
  end
  if complexity
    return (phenos,avg_complexity)
  else
    return phenos
  end
end
