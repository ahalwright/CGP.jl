# Uses map-reduce to union the set results of multiple runs of pheno_set_rand_neutral_walk() for a single phenotype ph.
# unions of the results of multiple walks
function pheno_rand_neutral_walks( p::Parameters, funcs::Vector{Func}, ph::Goal, nwalks::Int64, walk_length::Int64,  max_tries::Int64, max_steps::Int64;
    output_vect::Bool=false, use_lincircuit::Bool=false, csvfile::String="" )
  # One list of common phenotypes included here for reference.
  #common_list = [ 0x0000, 0x0011, 0x0022, 0x0033, 0x0044, 0x0055, 0x0066, 0x0077, 0x0088, 0x0099, 0x00aa, 0x00bb, 0x00cc, 0x00dd, 0x00ee, 0x00ff ]
  df = DataFrame()
  if output_vect
    pheno_vects  = pmap( _->pheno_rand_neutral_walk( p, funcs, ph, walk_length, max_tries, max_steps,
        output_vect=output_vect, use_lincircuit=use_lincircuit ), collect(1:nwalks ) ) 
    pheno_vect = reduce( +, pheno_vects )
    df.pheno_vect = pheno_vect
  else
    pheno_sets = pmap( _->pheno_rand_neutral_walk( p, funcs, ph, walk_length, max_tries, max_steps, 
        output_vect=output_vect, use_lincircuit=use_lincircuit ), collect(1:nwalks ) )
    pheno_set = reduce( union, pheno_sets )
    df.pheno_list =  map( x->@sprintf("0x%04x",x), [ s for s in pheno_set ] )
  end
  if length(csvfile) > 0
    hostname = chomp(open("/etc/hostname") do f read(f,String) end)
    open( csvfile, "w" ) do f
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      #println(f,"# run time in minutes: ",(ptime+ntime)/60)
      print_parameters(f,p,comment=true)
      println(f,"# funcs: ", Main.CGP.default_funcs(p))
      println(f,"# phenotype: ", ph)
      println(f,"# nwalks: ",nwalks)
      println(f,"# walk_length: ",walk_length)
      println(f,"# max_tries: ",max_tries)
      println(f,"# max_steps: ",max_steps)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  df
end

# returns the set or occurence vector of phenotypes encountered in a random neutral walk starting from a circuit evolved to map to phenotype ph.
function pheno_rand_neutral_walk( p::Parameters, funcs::Vector{Func}, ph::Goal, walk_length::Int64,  max_tries::Int64, max_steps::Int64;
      output_vect::Bool=false, use_lincircuit::Bool=false )
  (nc,steps) = pheno_evolve( p, funcs, ph::Goal, max_tries, max_steps; use_lincircuit=use_lincircuit ) 
  if steps == max_steps
    println("pheno_evolve failed to evolve a circuit to map to phenotype: ", ph )
    return Set(MyInt[])
  end
  pheno_rand_neutral_walk( nc, funcs, walk_length, max_tries, max_steps, output_vect=output_vect )
end


# returns the set or occurence vector of phenotypes encountered in a random neutral walk starting from a circuit evolved to map to phenotype ph.
function pheno_rand_neutral_walk( c::Circuit, funcs::Vector{Func}, walk_length::Int64,  max_tries::Int64, max_steps::Int64;
      output_vect::Bool=false )
  use_lincircuit = typeof(c) == LinCircuit
  println("use_lincircuit: ",use_lincircuit)
  p = c.params
  @assert p.numoutputs == 1
  goal = output_values(c)[1]
  if output_vect
    pheno_vect = zeros(Int64,2^2^p.numinputs)  # Assumes 1 output
  else
    pheno_set = Set(MyInt[])  # Assumes 1 output
  end
  for i = 1:walk_length
    j = 1
    while j <= max_steps   # also terminated by a break statement
      sav_c = deepcopy(c)
      use_lincircuit ? mutate_circuit!( c, funcs ) : mutate_chromosome!( c, funcs )
      #println("c: ",c)
      output = output_values(c)[1]
      if output == goal
        #println("successful step for i= ",i,"  outputs: ",outputs)
        all_neighbor_phenos = mutate_all( c, funcs )   # polymorphic call
        #println("j: ",j,"  all_neighbor_phenos: ",all_neighbor_phenos)
        pheno_list = map( x->x[1], all_neighbor_phenos )
        #println("j: ",j,"  pheno_list: ",pheno_list)
        if output_vect
          map( v->pheno_vect[v+1] += 1, pheno_list )
        else
          union!(pheno_set, Set(pheno_list) )
        end
        break
      end
      c = sav_c
      j += 1
    end
    if step == max_steps
      error("Failed to find neutral mutation")
    end
  end
  if output_vect
    pheno_vect
  else
    pheno_set
  end
end
#=
# returns the set of phenotypes encountered in a random neutral walk starting from a circuit evolved to map to phenotype ph.
function pheno_set_rand_neutral_walk( c::Circuit, funcs::Vector{Func}, walk_length::Int64,  max_tries::Int64, max_steps::Int64 )
  use_lincircuit = typeof(c) == LinCircuit
  println("use_lincircuit: ",use_lincircuit)
  p = c.params
  @assert p.numoutputs == 1
  goal = output_values(c)[1]
  pheno_set = Set(MyInt[])  # Assumes 1 output
  for i = 1:walk_length
    j = 1
    while j <= max_steps   # also terminated by a break statement
      sav_c = deepcopy(c)
      use_lincircuit ? mutate_circuit!( c, funcs ) : mutate_chromosome!( c, funcs )
      #println("c: ",c)
      output = output_values(c)[1]
      if output == goal
        #println("successful step for i= ",i,"  outputs: ",outputs)
        all_neighbor_phenos = mutate_all( c, funcs )   # polymorphic call
        #println("j: ",j,"  all_neighbor_phenos: ",all_neighbor_phenos)
        pheno_list = map( x->x[1], all_neighbor_phenos )
        #println("j: ",j,"  pheno_list: ",pheno_list)
        union!(pheno_set, Set(pheno_list) )
        break
      end
      c = sav_c
      j += 1
    end
    if step == max_steps
      error("Failed to find neutral mutation")
    end
  end
  pheno_set
end
=#
