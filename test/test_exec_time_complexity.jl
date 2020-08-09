# Measure execution time of a complexity, degeneracy, or redundancy function
# To run:  @time exec_time(10, complexity5, 4,4,16)  where complexity5 can be replaced by degeneracy or redundancy

maxints_for_degen = 20   #  Default 
function exec_time( reps::Int64, cmplx_fcn::Function, numinputs::Int64, numoutputs::Int64, numinteriors::Int64 )
  p = Parameters(numinputs=numinputs, numoutputs=numinputs, numinteriors=numinteriors, numlevelsback=numinteriors)
  funcs=default_funcs(numinputs)
  for _  = 1:reps
    c = random_chromosome( p, funcs )
    cc = cmplx_fcn( c )
  end
end

#  julia> @time exec_time(10, complexity5, 4,4,16)
#    1.979641 seconds (31.73 M allocations: 1.084 GiB, 5.20% gc time)
#  
#  julia> @time exec_time(10, degeneracy, 4,4,16)
#    4.620361 seconds (86.13 M allocations: 2.730 GiB, 4.38% gc time)
# Conclusion:  complexity5 is much faster than complexity4, complexity7
# Conclusion:  complexity5 is somewhat faster than complexity6 
# Conclusion:  degeneracy is much faster than degeneracy1

