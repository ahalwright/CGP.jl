# Note requires function mutational_robustness() which is currently in file Avg_mut_robustness
#(max_maxfit_robust, maxfit_chroms_list, max_ch_match, ch_match_indices, new_mut_robustness ) = Next_chromosome!(c,goallist,funcs,0.0,hamming_sel=true)
#(new_c,max_ch_fitness,max_robustness,goals_matched_list) = Next_chromosome!(c,goallist,funcs,0.0,hamming_sel=true)
export Next_chromosome!

function Next_chromosome!(c::Chromosome, goallist::GoalList, funcs::Vector{Func};
      hamming_sel::Bool=true, active_only::Bool=false, robust_select::Bool=true  )
  context = construct_context(c.params.numinputs)
  goals_matched = hamming_sel ? goals_matched_hamming : goals_matched_exact
  new_c = Chromosome[]
  new_component_score = zeros(Float64,c.params.lambda)
  maxfit_chroms_list = Tuple[]
  #maxfit_chroms_list = Float64[]
  max_new_component_score_ind = 1
  max_new_mut_robustness_ind  = 1
  for i = 1:c.params.lambda
    #print("i: ",i,"  ")
    push!(new_c, deepcopy(c))
    mutate_chromosome!( new_c[end], funcs )
    output = execute_chromosome(new_c[end],context)
    (best_score,best_score_list) = goals_matched( output, goallist, c.params.numinputs)
    #println("i: ",i,"  (best_score,best_score_list): ",(best_score,best_score_list))
    push!(maxfit_chroms_list,(best_score,best_score_list))
    #push!(maxfit_chroms_list,best_score)
  end
  #println("maxfit_chroms_list: ",maxfit_chroms_list)
  (max_ch_fitness, max_ch_fitness_indices)  = findmaxall_1( maxfit_chroms_list )
  #(max_ch_fitness, max_ch_fitness_indices)  = findmaxall( maxfit_chroms_list )
  #println("(max_ch_fitness, max_ch_fitness_indices):",(max_ch_fitness, max_ch_fitness_indices))
  if !robust_select
    max_ch_fitness_index = rand(max_ch_fitness_indices)  # randomly select the index for max_ch_fitness 
    return ( new_c[max_ch_fitness_index], max_ch_fitness[1], 0.0, maxfit_chroms_list[max_ch_fitness_index][2] )
  end  
  new_mut_robustness = Float64[]
  for j = 1:length(max_ch_fitness_indices)
    mut_robustness = mutational_robustness( c, funcs, active_only=active_only )
    #println("j: ",j,"  mut_robustness: ",mut_robustness)
    push!(new_mut_robustness,mut_robustness)
  end
  #println("new_mut_robustness: ",new_mut_robustness)
  (max_maxfit_robust, max_maxfit_robust_indices) = findmaxall( new_mut_robustness )  # best robustness out of max fitness chromosomes
  #println("max_maxfit_robust: ",max_maxfit_robust,"  max_maxfit_robust_indices: ",max_maxfit_robust_indices)  
  max_robust_indices = max_ch_fitness_indices[max_maxfit_robust_indices] # best robustness indices out of range 1:lambda 
  #println("max_robust_indices: ",max_robust_indices)# best robustness out of range 1:lambda
  max_chrome_index =   rand(max_robust_indices)  # randomly select a max robustness index
  #println("max_chrome_index: ",max_chrome_index,"  maxfit_chroms_list[max_chrome_index]: ",maxfit_chroms_list[max_chrome_index])
  max_robustness = maximum(new_mut_robustness)
  #println("max_robustness: ",max_robustness)
  goals_matched_list = maxfit_chroms_list[max_chrome_index][2]
  #println("goals_matched_list: ",goals_matched_list)
  #(max_maxfit_robust, maxfit_chroms_list, max_ch_fitness, max_ch_fitness_indices, new_mut_robustness )
  (new_c[max_chrome_index],max_ch_fitness[1],max_robustness,goals_matched_list)
end

# argument A is a vector where > and == can be used to compare elements
function findmaxall( A::AbstractVector )
  if isempty(A)
    error("Empty list in findmaxall()")
  end
  indices = [1]
  for i = 2:length(A) 
    if A[i] > A[indices[1]]
      indices = [i]
    elseif A[i] == A[indices[1]]
      push!(indices,i)
    end
  end
  (A[indices[1]],indices)
end
 
# argument A is a vector where < and == can be used to compare elements
function findminall( A::AbstractVector )
  if isempty(A)
    error("Empty list in findminall()")
  end
  indices = [1]
  for i = 2:length(A) 
    if A[i] < A[indices[1]]
      indices = [i]
    elseif A[i] == A[indices[1]]
      push!(indices,i)
    end
  end
  (A[indices[1]],indices)
end

# cmp is a two-argument comparison function:  greater-than for max, less-than for min
# eq  is a two-argument equality function
function findmaxminall( cmp, eq, A::AbstractVector )
  if isempty(A)
    error("Empty list in findmaxminall()")
  end
  indices = [1]
  for i = 2:length(A) 
    if cmp(A[i], A[indices[1]])
      indices = [i]
    elseif eq(A[i], A[indices[1]])
      push!(indices,i)
    end
  end
  (A[indices[1]],indices)
end

# findmaxall where the comparison functions look at the first components of the elements of A
findmaxall_1( A::AbstractVector ) = findmaxminall( (x,y)->(x[1]>y[1]), (x,y)->(x[1]==y[1]), A )

# findmnall where the comparison functions look at the first components of the elements of A
findminall_1( A::AbstractVector ) = findmaxminall( (x,y)->(x[1]<y[1]), (x,y)->(x[1]==y[1]), A )
 
