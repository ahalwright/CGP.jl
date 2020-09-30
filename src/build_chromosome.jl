

# Build the chromosome that corresponds to genes
function build_chromosome( p::Parameters, genes::Vector{Int64}, funcs::Vector{Func} )
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
      @assert genes[gi]+1 <= p.numlevelsback
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
    
# assumes that all funcs used have arity p.nodearity
function build_all_genes( p::Parameters )
  num_genes = p.numinteriors*p.nodearity
  max_gene = [ min(p.numlevelsback-1, i+p.numinputs-2) for i = 1:p.numinteriors for j=1:p.nodearity]
  println("max_gene: ",max_gene)
  current_genes = fill(0,num_genes)
  genes_list = [deepcopy(current_genes)]
  ints_list = [genes_to_integer(p,current_genes)]
  count = 0
  gi = 1
  while gi < num_genes 
    if current_genes[gi] < max_gene[gi]
      current_genes[gi] += 1
      #println("gi: ",gi,"  current_genes: ",current_genes)
      push!(genes_list,deepcopy(current_genes))
      push!(ints_list,genes_to_integer(p,current_genes))
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
        push!(genes_list,deepcopy(current_genes))
        push!(ints_list,genes_to_integer(p,current_genes))
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
  (genes_list, ints_list)
end

function chromosome_to_genes( c::Chromosome )
  genes = zeros(Int64,p.nodearity*p.numinteriors)
  ints = c.interiors
  for i = 1:length(ints)
    for j = p.nodearity:-1:1
      genes[ p.nodearity*(i-1) + j] = i-ints[i].inputs[j]+1
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
    
function build_outputs_list( p::Parameters, funcs::Vector{Func}, genes_list::Vector{Vector{Int64}} )
  result = fill(MyInt[],length( genes_list ))
  for g in genes_list
    index = genes_to_integer( p, g )
    c=build_chromosome(p,g,funcs)
    output = output_values( c )
    #print_build_chromosome(c)
    result[index+1] = output
  end
  result
end  
