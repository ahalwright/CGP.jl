# Construct the phenotype network for a given parameter setting.
using HDF5, JLD

function create_empty_pheno_net( p::Parameters; pheno_file::String="" )
  pheno_net = zeros( Int64, p.numoutputs*2^2^p.numinputs, p.numoutputs*2^2^p.numinputs )
  if pheno_file != ""
    save( pheno_file, "pheno_net", pheno_net )
  end
  pheno_net
end

function construct_pheno_net( p::Parameters, nreps::Int64, numcircuits::Int64; pheno_file::String="", use_lincircuit::Bool=false )
  funcs = default_funcs(p.numinputs)
  if pheno_file == "" || !isfile(pheno_file)
    ph_net = create_empty_pheno_net( p )
  elseif isfile(pheno_file)
    ph_net = jldopen(pheno_file,"r") do file
      read( file, "pheno_net" )
    end 
  end
  (outlist, ch_ints_list) = count_outputs( nreps, p, numcircuits, use_lincircuit=use_lincircuit )
  println("length(outlist): ",length(outlist))
  println("length(ch_ints_list): ",length(ch_ints_list))
  for ch_ints in ch_ints_list
    println("ch_ints_list: ",ch_ints_list)
    for ch_int in ch_ints
      println("ch_int: ",ch_int)
      if use_lincircuit
        ch_vects = instruction_ints_to_instruction_vects( circuit_int_to_instruction_ints( ch_int, p, funcs ), p, funcs )  
        ch = LinCircuit( ch_vects, p ) 
      else
        #ch = ch_int.circuit_vects
        ch = int_to_chromosome( ch_int, p, funcs )
        #print_circuit(ch)
      end
      src = output_values(ch)[1]   # assumes 1 output    
      println("ch_int: ",ch_int,"  src: ",src)
      dests = map(x->x[1], mutate_all( ch, funcs ))
      println("dests: ",dests)
      for dest in dests
        ph_net[src+1,dest+1] += 1
        #println("dest: ",dest,"  ph_net[src+1,dest+1]: ",ph_net[src+1,dest+1])
      end
    end
  end
  if pheno_file != ""
    save( pheno_file, "pheno_net", ph_net )
  end
  ph_net
end
