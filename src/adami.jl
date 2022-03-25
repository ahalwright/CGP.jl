# Circuit entropy and Adami complexity

# based on random_chromosome() in Chromosome.jl
# The three components of an interior node of a chromosome are the function and the two inputs.
# Each of these is specified by an integer.  
# The corresponding arity is the pair of the number of possible values for the integer and the minimum value of the integer.
# Thus, the arity of an interior node is the triple of these pairs.
# gate_arities() returns the list of arities for each of the interior nodes of a chromosome corresponding to paramters p.
function gate_arities( p::Parameters, funcs::Vector{Func} )
  gate_arities = Vector{Tuple{Int64,Int64}}[]
  for index = (p.numinputs +1):(p.numinputs+p.numinteriors)
    maxindex = index - 1
    minindex = max(1,index-p.numlevelsback) 
    arities = [(length(funcs),1),(maxindex-minindex+1,minindex-1),(maxindex-minindex+1,minindex-1)]
    push!(gate_arities,arities)
  end
  gate_arities
end

# A record of the number of chromosomes that have a value for each of these fields
mutable struct GateValues
  func::Vector{Int64}
  input1::Vector{Int64}
  input2::Vector{Int64}
end

# Returns a list of zero GateValues, one for each of the p.numinteriors interior nodes.
# Each GateValues vector is a zeros vector whose length is specified by arities[1] for that interior node.
function init_chrome_values( p::Parameters, funcs::Vector{Func} )
  arities = gate_arities( p, funcs )
  chrome_values = GateValues[]
  for i = 1:p.numinteriors
    chv = GateValues( zeros(Int64,arities[i][1][1]), zeros(Int64,arities[i][2][1]), zeros(Int64,arities[i][3][1]))
    push!(chrome_values,chv)
  end
  chrome_values
end

# update chrome_values for one chromosome by incrementing one component of the corresponding chrome_values array by 1.
function update_chrome_values!( chrome_values::Vector{GateValues}, arities::Vector{Vector{Tuple{Int64, Int64}}}, ch::Chromosome, funcs::Vector{Func} )
  @assert ch.params.numinteriors <= length(chrome_values)
  #print_circuit(ch)
  for i = 1:ch.params.numinteriors
    j = 1  # determine index of func of interior node i
    while j <= length(funcs) && ch.interiors[i].func.func != funcs[j].func
      j += 1
    end
    #println("i: ",i,"  j: ",j)
    chrome_values[i].func[j] += 1
    chrome_values[i].input1[ch.interiors[i].inputs[1]-arities[i][2][2]] += 1
    chrome_values[i].input2[ch.interiors[i].inputs[2]-arities[i][3][2]] += 1
  end
end

# update chrome_values for a population of chromosomes
function update_chrome_values!( chrome_values::Vector{GateValues}, arities::Vector{Vector{Tuple{Int64, Int64}}}, pop::Vector{Chromosome}, funcs::Vector{Func} )
  for ch in pop
    update_chrome_values!( chrome_values, arities, ch, funcs )
  end
end

# The entropy of a vector of Integers where v[i] is the number of chromosomes with value i for an interior node field
# This is just the standard entropy formula
function aentropy( v::Vector{Int64} )
  log2z(x) = (x == 0.0 ? 0.0 : log2(x))  # log2 with 0 mapped to 0
  sumv = sum(v)
  if sumv != 0 
    p = v/sumv
  else
    return 0.0
  end
  sum( -p[i]*log2z(p[i]) for i = 1:length(p) )
end

function entropy( gv::Vector{GateValues} )
  entropies = Float64[]
  for g in gv
    ents = [aentropy(g.func),aentropy(g.input1),aentropy(g.input2)]
    entropies = vcat(entropies,ents)
  end
  entropies
end
  
function pop_entropies( p::Parameters, funcs::Vector{Func}, pop::Vector{Chromosome} )
  for ch in pop
    @assert ch.params.numinputs == p.numinputs
    @assert ch.params.numoutputs == p.numoutputs
    @assert ch.params.numinteriors == p.numinteriors
    @assert ch.params.numlevelsback == p.numlevelsback
  end
  arities=gate_arities(p,funcs)
  chvalues = init_chrome_values( p, funcs )
  update_chrome_values!( chvalues, arities, pop, funcs )
  #println("chvalues: ",chvalues)
  entropy(chvalues)
end

function avg_pop_entropy( pop::Vector{Chromosome}, funcs::Vector{Func}=default_funcs(pop[1].params) )
  arities=gate_arities(p,funcs)
  chvalues = init_chrome_values( p, funcs )
  update_chrome_values!( chvalues, arities, pop, funcs )
  #println("chvalues: ",chvalues)
  entvalues = entropy(chvalues)
  sum(entvalues)/length(entvalues)
end

# Tests only that there are no indexing errors
function test_update_chrome_values(p,popsize::Int64,funcs=default_funcs(p))
  pop = Chromosome[]
  for i = 1:popsize
    push!(pop,random_chromosome(p,funcs))
  end
  arities=gate_arities(p,funcs)
  chvalues = init_chrome_values( p, funcs )
  update_chrome_values!( chvalues, arities, pop, funcs )
  chvalues
end

#=
p = Parameters(2,1,4,3)
funcs = default_funcs(p)
pop = Chromosome[]
println("arities: ",circuit_arities(p,funcs))
popsize = 4
for i = 1:popsize
  push!(pop,random_chromosome(p,funcs))
end
circuit_pop_entropy( p, funcs, pop )
for j = 1:popsize
  println(ch_values(pop[j],3 ))
end
=#
