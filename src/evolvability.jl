# Define phenotypic evolvability where a phenotype is a goal or chromosome output.
# Note that genotypic evolvability is computed by calling mutate_all( c, funcs, robust_only=true )[2]
# Use Wagner's (2008) method of evolving nchromss chromsomes whose output is the goal (where Wagner uses nchroms=100).

function evolvability( g::Goal, funcs, nchroms,  maxsteps, numinputs, numinteriors, numlevelsback )
  p = Parameters( numinputs=numinputs, numoutputs=length(goal), numinteriors=numinteriors, numlevelsback=numlevelsback )
  result = Goal[]
  for i = 1:nchroms
    c = random_chromosome( p, funcs )
    c = mut_evolve( c, [g], funcs, maxsteps )[1]
    print_build_chromosome( c )
    goal_list = mutate_all( c, funcs, output_outputs=true )
    #println("goal_list: ",goal_list)
    goal_set = unique(goal_list)
    println("goal_set: ",goal_set)
    result = unique(vcat(result,goal_set))
  end
  result
end
    
# mut_evolve( c::Chromosome, goallist::GoalList, funcs::Vector{Func}, max_steps::Integer;    
# function mutate_all( c::Chromosome, funcs::Vector{Func};
#      robustness_only::Bool=false, output_outputs::Bool=true, output_chromosomes::Bool=false )
    
