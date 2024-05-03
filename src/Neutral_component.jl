using DataFrames
# returns the component of the neutral set of the phenotype output_values(circuit) that contains circuit.
# If the evolvability keyword argument is true, also returns a list of the phenotypes seen.
# Tested by the function test_neutral_components() in test/test_neutral_componemts.jl
# Example:  p = Parameters(2,1,3,3); funcs=default_funcs(p)  # 5 funcs
# ch = int_to_chromosome( 39403, p, funcs )
# nc0 = neutral_component( ch, funcs )
# Set{Int128} with 4332 elements
# Consistent with the results of component_properties() from Fnc_mt.jl on 4/28/24
function neutral_component( ch::Chromosome, funcs::Vector{Func}, evolvability::Bool=false )::Union{Set,Tuple{Set,Vector{Goal}}}
  p = ch.params
  ov = output_values(ch)
  nc = Set([chromosome_to_int(ch,funcs)])
  phenos = Goal[ov]  
  stack = [chromosome_to_int(ch,funcs)]
  while length(stack) > 0  #&& cnt > 0
    circ = int_to_chromosome( pop!(stack), p, funcs )
    circs = filter(cch->ov==output_values(cch), mutate_all( circ, funcs, output_outputs=false, output_circuits=true ))
    #println("length(stack): ",length(stack),"  length(circs): ",length(circs))
    for c in circs
      if !( chromosome_to_int(c) in nc )
        push!(nc,chromosome_to_int(c,funcs))
        push!(stack,chromosome_to_int(c,funcs))
        if evolvability
          append!(phenos,mutate_all(c,funcs,output_outputs=true,output_circuits=false))
        end
      end
    end
  end
  if evolvability
    return (nc,unique(phenos))
  else
    return nc
  end
end

function neutral_component( ch::LinCircuit, funcs::Vector{Func}, evolvability::Bool=false )
  p = ch.params
  ov = output_values(ch)
  nc = Set([circuit_to_circuit_int(ch,funcs)])
  phenos = Goal[ov]  
  stack = [circuit_to_circuit_int(ch,funcs)]
  cnt = 70
  while length(stack) > 0 && cnt > 0
    #println("length(stack): ",length(stack))
    circ = circuit_int_to_circuit( pop!(stack), p, funcs )
    circs = filter(cch->ov==output_values(cch), mutate_all( circ, funcs, output_outputs=false, output_circuits=true ))
    for c in circs
      if !( circuit_to_circuit_int(c,funcs) in nc )
        push!(nc,circuit_to_circuit_int(c,funcs))
        push!(stack,circuit_to_circuit_int(c,funcs))
        if evolvability
          append!(phenos,mutate_all(c,funcs,output_outputs=true,output_circuits=false))
        end
      end
    end
    cnt -= 1
  end
  if evolvability
    return (nc,unique(phenos))
  else
    return nc
  end
end

# Evolves a chromosome for each phenotype and then computes the lengths of the neutral components for each of these chromosomes.
function count_neutral_components_ch_mt( p::Parameters, funcs::Vector{Func} )
  ncomponent_counts = [ Threads.Atomic{Int64}(0) for i= 1:2^2^p.numinputs ]
  Threads.@threads for ph = MyInt(0):MyInt(2^2^p.numinputs-1)
    rch = random_chromosome( p, funcs )
    ch_ph = neutral_evolution( rch, funcs, [ph], 30_000 )[1]
    #print_circuit( ch_ph )
    #println("output_values: ",output_values(ch_ph))
    ncomponents = neutral_component( ch_ph, funcs )
    #println("ph: ",ph,"  ncomponents: ",length(ncomponents))
    Threads.atomic_add!( ncomponent_counts[ph+1], length(ncomponents) )
  end
  gc = map( ph->ph[], ncomponent_counts )
end

# Evolves a chromosome for each phenotype and then computes the lengths of the neutral components for each of these chromosomes.
function count_neutral_components_ch( p::Parameters, funcs::Vector{Func} )
  #ncomponent_counts = [ Threads.Atomic{Int64}(0) for i= 1:2^2^p.numinputs ]
  ncomponent_counts = zeros(Int64,2^2^p.numinputs)
  for ph = MyInt(0):MyInt(2^2^p.numinputs-1)
    rch = random_chromosome( p, funcs )
    ch_ph = neutral_evolution( rch, funcs, [ph], 30_000 )[1]
    #print_circuit( ch_ph )
    #println("output_values: ",output_values(ch_ph))
    ncomponents = length(neutral_component( ch_ph, funcs ))
    #println("ch_ph: ",ch_ph,"  ncomponents: ",ncomponents)
    ncomponent_counts[ph+1] = ncomponents
  end
  ncomponent_counts
end

# returns a DataFrame with one row per phenotype.  The :count column is the number of genotypes for each phenotype.
# Then there are n_columns additional columns which are the results of function count_neutral_components_ch_mt() defined above
function neutral_components_dataframe_mt( p::Parameters, funcs::Vector{Func}, n_columns::Int64=1 )::DataFrame
  df = count_genos = count_genotypes_ch_mt( p, funcs )
  #count_ncomponents = count_neutral_components_ch( p, funcs )
  for i = n_columns:-1:1
    column_symbol = Symbol("ncomponents$(i)")
    insertcols!( df, 2, column_symbol=>count_neutral_components_ch_mt( p, funcs ) )
  end
  df
end

#function neutral_component_walk_dataframe_mt( p::Parameters, funcs::Vector{Func}, n_columns::Int64=1 )::DataFrame

function neutral_components_df_to_csv( p::Parameters, funcs::Vector{Func}, df::DataFrame, csvfile::String="" )
  if length(csvfile) > 0
    println("csvfile: ",csvfile)
    open( csvfile, "w" ) do f
      hostname = readchomp(`hostname`)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# funcs: ", funcs )
      print_parameters(f,p,comment=true)
      #println(f,"# n_columns: ",n_columns)
      CSV.write( f, df, append=true, writeheader=true )
      println("csvfile written")
    end
  end
end

# The set components is a set of Chromosomes.
# If this function doesn't stack overflow, it returns the neutral component of the neutral set of circuit.
# The Chromsome circuit is added to components along with all mutations of circuit.
# If this did not increase the size of components, quit.
# Otherwise, recursively call this function on all of these mutations.
# The above function is both more efficient and more reliable:   outdated
function neutral_component_ints!( components::Set,  circuit::Chromosome )
  size0 = length( components )
  #println("size0: ",size0)
  ov = output_values(circuit)
  mcn = filter(ch->ov==output_values(ch), mutate_all( circuit, funcs, output_outputs=false, output_circuits=true ))
  mcn_ints = map( c->chromosome_to_int( c, funcs ), mcn )
  union!( components, Set( mcn_ints ) )
  new_size = length( components )
  if size0 == new_size
    return components
  else
    for c in mcn
      neutral_component_ints!( components, c )
    end
  end
end

# The set components is a set of Chromosomes.
# If this function doesn't stack overflow, it returns the neutral component of the neutral set of circuit.
# The Chromsome circuit is added to components along with all mutations of circuit.
# If this did not increase the size of components, quit.
# Otherwise, recursively call this function on all of these mutations.
function neutral_component_ints( circuit::Chromosome )
  components = Set([])
  neutral_component_ints!( components::Set,  circuit )
end
