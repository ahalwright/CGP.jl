using DataFrames, CSV, Printf, Distributed
export build_all_genes, build_chromosome, chromosome_to_genes, genes_to_integer, integer_to_genes
export build_outputs_list, outputs_list,  build_outputs, save_ints_goals_to_file, build_max_gene
    
# assumes that all funcs used have arity p.nodearity
#function build_all_genes( p::Parameters; max_gene::Vector{Int64}=Int64[] )
function build_all_genes( p::Parameters; incr::Int64=0 )
  num_genes = p.numinteriors*p.nodearity
  max_gene = build_max_gene(p,incr)
  #=
  if length(max_gene) == 0
    max_gene = [ min(p.numlevelsback-1, i+p.numinputs-2) for i = 1:p.numinteriors for j=1:p.nodearity]
  end
  =#
  println("max_gene: ",max_gene)
  current_genes = fill(0,num_genes)
  pairs_list = [ (genes_to_integer(p,current_genes), deepcopy(current_genes)) ]
  #genes_list = [deepcopy(current_genes)]
  #ints_list = [genes_to_integer(p,current_genes)]
  count = 0
  gi = 1
  while gi < num_genes 
    if current_genes[gi] < max_gene[gi]
      current_genes[gi] += 1
      #println("gi: ",gi,"  current_genes: ",current_genes)
      push!(pairs_list, (genes_to_integer(p,current_genes),deepcopy(current_genes)))   
      #push!(genes_list,deepcopy(current_genes))
      #push!(ints_list,genes_to_integer(p,current_genes))
      #println("glen: ",length(genes_list),"  ilen: ",length(ints_list))
    else
      gi += 1
      while gi < num_genes && current_genes[gi] == max_gene[gi]
        gi += 1
      end
      if current_genes[gi] < max_gene[gi]
        current_genes[gi] += 1
        current_genes[1:(gi-1)] = fill(0,(gi-1))
        #println("gi: ",gi,"  current_genes: ",current_genes)
        push!(pairs_list, (genes_to_integer(p,current_genes),deepcopy(current_genes)))   
        #push!(genes_list,deepcopy(current_genes))
        #push!(ints_list,genes_to_integer(p,current_genes))
      else
        break
      end
      gi = 1
    end
    #println("gi: ",gi,"  cg: ",current_genes)
    count += 1
    if count % 10000000 == 0
      println("count: ",count)
    end
  end
  #(genes_list, ints_list)
  pairs_list
end

function build_max_gene( p::Parameters, incr::Int64 )
  max_gene = [ min(p.numlevelsback-1, i+p.numinputs-2) for i = 1:p.numinteriors for j=1:p.nodearity]
  for j=1:p.nodearity
    for i = 1:p.numinteriors
      if (j-1)%p.nodearity == 0
        max_gene[(2*i-1)+(j-1)] = min(p.numlevelsback+incr-1, i+p.numinputs-2)
      end
    end
  end
  max_gene
end

# Build the chromosome that corresponds to genes
function build_chromosome( p::Parameters, genes::Vector{Int64}, funcs::Vector{Func}; incr::Int64=0 )
  @assert p.nodearity == 2
  inputs = Vector{InputNode}[]
  interiors = Vector{InteriorNode}[]
  outputs = Vector{OutputNode}[]
  c = Chromosome(p)
  for index = 1:p.numinputs
    c.inputs[index] = InputNode(index)
    c.inputs[index].active = false
    c.inputs[index].cache = 0
  end
  gene_index = 1
  for index = (p.numinputs +1):(p.numinputs+p.numinteriors)
    func = funcs[rand(1:end)]  # replace later
    @assert func.arity == 2
    inputs = zeros(Int64,func.arity)
    i = 1
    for gi = gene_index:(gene_index+func.arity-1)
      #print("index: ",index,"  gi: ",gi,"  genes[gi]: ",genes[gi],"  i: ",i)
      @assert genes[gi]+1 <= p.numlevelsback+incr
      @assert index-(genes[gi]+1) >= 1
      inputs[i] = index-(genes[gi]+1)
      #println(" inputs[i]: ",inputs[i])
      i += 1
    end
    #inputs = [ index-genes[gi] for gi = gene_index:(gene_index+func.arity-1) ]
    #println("index: ",index,"  gene_index: ",gene_index,"  inputs: ",inputs)
    gene_index += func.arity
    c.interiors[index - p.numinputs] = InteriorNode(func, inputs)
  end
  for i = 1:p.numoutputs
    index = p.numinputs + p.numinteriors + i - p.numoutputs   # use the last numoutputs interiors
    c.outputs[i] = OutputNode(index)
    c[index].active = true  #  Output nodes are always active
  end
  set_active_to_false(c)
  return c
end

function chromosome_to_genes( c::Chromosome )
  genes = zeros(Int64,c.params.nodearity*c.params.numinteriors)
  ints = c.interiors
  for i = 1:length(ints)
    for j = c.params.nodearity:-1:1
      genes[ c.params.nodearity*(i-1) + j] = i-ints[i].inputs[j]+1
    end
  end
  genes
end

# Note that result can be zero, so when used as an array index, add 1
function genes_to_integer( p::Parameters, genes::Vector{Int64} )
  num_genes = p.numinteriors*p.nodearity
  max_gene = [ min(p.numlevelsback-1, i+p.numinputs-2) for i = 1:p.numinteriors for j=1:p.nodearity]
  #println("max_gene: ",max_gene)
  result = 0
  multiplier = 1
  for i = num_genes:-1:1
    result += multiplier*genes[i]
    #println("i: ",i,"  result: ",result)
    multiplier *= max_gene[i]+1
    #println("i: ",i,"  multiplier: ",multiplier)
  end
  result
end

function integer_to_genes( p::Parameters, intg::Int64 )
  num_genes = p.numinteriors*p.nodearity
  max_gene = [ min(p.numlevelsback-1, i+p.numinputs-2) for i = 1:p.numinteriors for j=1:p.nodearity]
  #println("max_gene: ",max_gene)
  genes = zeros(Int64,num_genes)
  ig = intg
  for i = num_genes:-1:1
    genes[i] = ig % (max_gene[i]+1)
    #println("i: ",i,"  genes[i]: ",genes[i])
    ig = Int( floor( ig/(max_gene[i]+1) ) )
    #println("ig: ",ig) 
  end
  genes
end
    
function build_outputs_list( p::Parameters, funcs::Vector{Func}, int_gene_pairs::Vector{Tuple{Int64,Vector{Int64}}} )
  result = fill(MyInt[],length( int_gene_pairs ))
  index = 1
  for ig in int_gene_pairs
    c=build_chromosome(p,ig[2],funcs)
    output = output_values( c )
    #print_build_chromosome(c)
    result[index] = output
    index += 1
  end
  result
end 

# parallel version which splits gene_int_list into pieces of length piece_length for parallelization
function outputs_list( p::Parameters, int_gene_pairs::Vector{Tuple{Int64,Vector{Int64}}}, funcs::Vector{Func}, piece_length::Int64;
    incr::Int64=0 )
  upper_bound = Int(floor(length(int_gene_pairs)/piece_length))
  len = length(int_gene_pairs)
  println("len(gil): ",length(int_gene_pairs),"  upper_bound: ",upper_bound)
  #println( [ (i*piece_length+1,(min(i*piece_length+piece_length,len))) for i = 0:upper_bound ] )
  gene_sublists_list = [ int_gene_pairs[(i*piece_length+1):min(i*piece_length+piece_length,len)] for i = 0:upper_bound ]
  #result = map( g->build_outputs( p, g, funcs, incr=incr ), gene_sublists_list )
  result = pmap( g->build_outputs( p, g, funcs, incr=incr ), gene_sublists_list )
  vcat( result... )   # combine all of the results together
end


function build_outputs( p::Parameters, gene_int::Int64, funcs::Vector{Func} )
  c=build_chromosome(p,integer_to_genes(p,gene_int),funcs)
  (gene_int, output_values( c ))
end

# 
function build_outputs( p::Parameters, int_gene_pairs::Vector{Tuple{Int64,Vector{Int64}}}, funcs::Vector{Func}; incr::Int64=0  )
  result = Vector{Tuple{Int64, Vector{MyInt}}}[]
  result = fill( (0,[MyInt(0x0)]),length(int_gene_pairs))
  for i = 1:length(int_gene_pairs)
    c=build_chromosome(p,int_gene_pairs[i][2],funcs,incr=incr)
    result[i] = (int_gene_pairs[i][1], output_values( c ))
  end
  result
end

# Saves as a dataframe
# Assumes chromosomes have only one output 
function save_ints_goals_to_file( p::Parameters, int_gene_pairs::Vector{Tuple{Int64,Vector{MyInt}}}, csvfile::String; incr::Int64=0 )
  df = DataFrame( ch_ints=map( ig->ig[1], int_gene_pairs ), goals=map( x->@sprintf("0x%x",x[1][1]), int_gene_pairs ) )
  open(csvfile,"w") do f
    print_parameters(f,p,comment=true)
    if incr > 0
      println("incr: ",incr)
    end
    CSV.write( f, df, append=true, writeheader=true )
  end
end

function save_build_goals_to_file( bag::Vector{Tuple{Int64,Vector{Int64}}}, filename::String )
  open( filename, "w" ) do f
    print_parameters(f,p,comment=true)
    println(f,"# incr: ",incr)
    for i = 1:length(bag)
      println(f,bag[i])
    end
  end
end
