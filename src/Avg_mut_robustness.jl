# Implement experiments as described in the document "Research update 3_18_29.docx"
export avg_mutational_robustness, run_avg_mutational_robustness
using DataFrames 
using CSV
using Statistics

#=
numinputs = 2
numoutputs = 2
nodearity = 2
numinteriors = 6
numlevelsback = 4
context = construct_contexts(numinputs)[numinputs]
p = Parameters(numinputs, numoutputs, nodearity, numinteriors, numlevelsback)
=#
#iterations = 20
#print_parameters(p)

#=
# Computes the fraction of mutations that do not change the output
# If active_only==true, then only mutations that change active nodes are considered
# If active_only==false, then all mutations are considered
# Note that results are partially random since the result of mutations at a specific location can be different.
function mutational_robustness( c::Chromosome, funcs::Vector{Func}; active_only::Bool=false ) 
  context = construct_context( c.params.numinputs )
  prev_out = execute_chromosome( c, context )  
  count_no_change = 0
  count_no_change = 0.0
  num_mut_locs = 0
  for j = 1:num_mutate_locations( c, funcs )
    cc = deepcopy(c)
    (cc,active) = mutate_chromosome!( cc, funcs, j )
    out_c = execute_chromosome( cc, context )
    #println("j: ",j,"  active: ",active)
    #print_chromosome( cc )
    if !active_only || active # consider all mutations if !active_only, only active mutations if active_only
      #out_c = execute_chromosome( cc, context )
      count_no_change += ((out_c == prev_out) ? 1 : 0)
      if active_only && active
        num_mut_locs += 1
      end
    end
    #println("j: ",j,"  active: ",active,"  no_change: ",(out_c == prev_out),"  count_no_change: ",count_no_change)
  end
  #println("num_mut_locs: ",num_mut_locs,"  num_mutate_locations(): ",num_mutate_locations(c,funcs),"  count_no_change: ",count_no_change)
  #@assert num_mut_locs == num_mutate_locations(c,funcs)
  if !active_only
    num_mut_locs = num_mutate_locations( c, funcs )
  end
  count_no_change/num_mut_locs
end
=#

function run_avg_mutational_robustness( iterations::Integer, numinputs::AbstractRange{Int64},
    numoutputs::AbstractRange{Int64}, numinteriors::AbstractRange{Int64}, numlevelsback::AbstractRange{Int64},
    funcs::Vector{Func}, csvfile::String; active_only::Bool=false )
  nodearity = 2
  df = DataFrame() 
  df.numinputs=Int64[]
  df.numoutputs=Int64[]
  df.numinteriors=Int64[]
  df.numlevelsback=Int64[]
  df.mut_robustness=Float64[]
  df.redundancy=Float64[]
  df.complexity=Float64[]
  df.degeneracy=Float64[]
  df.fract_active=Float64[]
  for num_inputs = numinputs
    for num_outputs = numoutputs
      for num_interiors = numinteriors
        for num_levelsback = numlevelsback
          p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, num_levelsback )
          #print_parameters( p )
          c = random_chromosome( p, funcs )
          (iterations,robustness,redundancy,complexity,degeneracy,fract_active) = avg_mutational_robustness( c, p, funcs, 
              iterations, active_only=true )
          #println("robustness: ",robustness,"  fract_active: ",fract_active)
          new_row = (num_inputs,num_outputs,num_interiors,num_levelsback,robustness,redundancy,complexity,degeneracy,fract_active)
          Base.push!( df, new_row )
        end
      end
    end
  end
  df
  open( csvfile, "w" ) do f
    println(f,"funcs: ", funcs )
    println(f,"nodearity: ",nodearity)
    println(f,"iterations: ",iterations)
    CSV.write( f, df, append=true, writeheader=true )
  end
  println(df)
end

# Currently uses fract_active only in correlation
function avg_mutational_robustness( c::Chromosome, p::Parameters, funcs::Vector{Func}, iterations::Integer;
    active_only::Bool=false )
  #println("active_only: ",active_only," in:",p.numinputs," out:",p.numoutputs," ints:",p.numinteriors)
  base = 2.0
  context = construct_context( p.numinputs )
  sum_redundancy = 0.0
  sum_complexity = 0.0
  sum_degeneracy = 0.0
  sum_mutational_robustness = 0.0
  sum_fract_active = 0.0
  sum_num_mut_locations = 0
  # Lists are for the computation of correlations
  fract_active_list = fill(0.0,iterations)
  mutational_robustness_list = fill(0.0,iterations)
  for i = 1:iterations
    c = random_chromosome( p, funcs )
    redund = redundancy( c, base=base )
    complexity = complexity5( c, base=base )
    degen = degeneracy( c, base=base )
    fract_active = fraction_active(c)
    mut_robustness = mutational_robustness( c, funcs, active_only=active_only )
    sum_mutational_robustness += mut_robustness
    sum_redundancy += redund
    sum_complexity += complexity
    sum_degeneracy += degen
    sum_fract_active += fract_active
    #sum_num_mut_locations += num_mut_locs
    #println("i: ",i,"  average_fract_active: ",average_fract_active)
    fract_active_list[i] = fract_active
    #println("i: ",i,"  mutational_robustness: ",mutational_robustness)
    mutational_robustness_list[i] = mut_robustness
  end
  correlation = cor( mutational_robustness_list, fract_active_list )
  #println("correlation: ",correlation)
  #println("total mutational_robustness: ",sum_mutational_robustness/iterations)
  (iterations,sum_mutational_robustness/iterations,sum_redundancy/iterations,
      sum_complexity/iterations,sum_degeneracy/iterations,sum_fract_active/iterations)

end
