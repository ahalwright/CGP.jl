# Implement experiments as described in the document "Research update 3_18_29.docx"
export avg_mutational_robustness, run_avg_mutational_robustness
using DataFrames 
using CSV
using Statistics


#= Parameters are set in SetParams.jl
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

function run_avg_mutational_robustness( iterations::Integer, numinputs::AbstractRange{Int64},
    numoutputs::AbstractRange{Int64}, numinteriors::AbstractRange{Int64}, numlevelsback::AbstractRange{Int64},
    funcs::Vector{Func}, csvfile::String; active_only::Bool=false )
  nodearity = 2
  df = DataFrame() 
  df.numinputs=Int64[]
  df.numoutputs=Int64[]
  df.numinteriors=Int64[]
  df.numlevelsback=Int64[]
  df.iterations=Int64[]
  df.mut_robustness=Float64[]
  df.frac_active=Float64[]
  df.correlation=Float64[]  
  for num_inputs = numinputs
    for num_outputs = numoutputs
      for num_interiors = numinteriors
        for num_levelsback = numlevelsback
          p = Parameters( num_inputs, num_outputs, nodearity, num_interiors, num_levelsback )
          #print_parameters( p )
          c = random_chromosome( p, funcs )
          (iterations,robustness,fract_active,correlation) = avg_mutational_robustness( c, p, funcs, 
              iterations, active_only=true )
          #println("robustness: ",robustness,"  fract_active: ",fract_active)
          new_row = (num_inputs,num_outputs,num_interiors,num_levelsback,iterations,robustness,fract_active,correlation)
          Base.push!( df, new_row )
        end
      end
    end
  end
  df
  CSV.write( csvfile, df )
  println(df)
end

function avg_mutational_robustness( c::Chromosome, p::Parameters, funcs::Vector{Func}, iterations::Integer;
    active_only::Bool=false )
  println("active_only: ",active_only," in:",p.numinputs," out:",p.numoutputs," ints:",p.numinteriors)
  sum_mutational_robustness = 0.0
  sum_avg_fract_active = 0.0
  sum_num_mut_locations = 0
  fract_active_list = fill(0.0,iterations)
  mutational_robustness_list = fill(0.0,iterations)
  for i = 1:iterations
    c = random_chromosome( p, funcs )
    #println("i: ",i,"  c: ")
    #print_chromosome( c )
    count_no_change = 0
    #num_mut_locs = num_mutate_locations( c, funcs )
    goal = execute_chromosome( c, context )
    sum_fract_active = 0.0
    num_mut_locs = 0
    for j = 1:num_mutate_locations( c, funcs )
      cc = deepcopy(c)
      (cc,active) = mutate_chromosome!( cc, funcs, j )
      #println("j: ",j,"  cc: ")
      #print_chromosome( cc )
      out_c = execute_chromosome( cc, context )
      if active_only && active
        sum_fract_active += fraction_active(cc)
        num_mut_locs += 1
      else
        sum_fract_active += fraction_active(cc)
      end
      #println("sum_fract_active: ",sum_fract_active)
      if out_c == goal
        count_no_change += 1
      end
    end
    if !active_only
      num_mut_locs = num_mutate_locations( c, funcs )
    end
    #println("num_mut_locs: ",num_mut_locs,"  num_mutate_locations(c,funcs): ",num_mutate_locations(c,funcs))
    sum_num_mut_locations += num_mut_locs
    average_fract_active = sum_fract_active/num_mut_locs
    #println("i: ",i,"  average_fract_active: ",average_fract_active)
    fract_active_list[i] = average_fract_active
    sum_avg_fract_active += average_fract_active
    mutational_robustness = count_no_change/num_mutate_locations( c, funcs )
    #println("i: ",i,"  mutational_robustness: ",mutational_robustness)
    mutational_robustness_list[i] = mutational_robustness
    sum_mutational_robustness += mutational_robustness
  end
  correlation = cor( mutational_robustness_list, fract_active_list )
  #println("correlation: ",correlation)
  #println("total mutational_robustness: ",sum_mutational_robustness/iterations)
  (iterations,sum_mutational_robustness/iterations,sum_avg_fract_active/iterations,correlation)
end
  
