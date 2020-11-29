# test of sampling hypotheses for number of genotypes that map to a given phenotype.
#=
Test the following hypothesis:

Given an urn with $n-k$ black balls and $k$ white balls.
Take repeated samples with the following procedure:
  * pick a ball from the urn
  * if the ball is black, return to the urn
  * if the ball is white, replace with a black (or blue) ball
Modeling the probability $p_i$ of a white ball on sample $i$ for $i = 0,1,2,\ldots$.
Hypothesis:  $p_i = k/n*(1-1/n)^i$. 
=#


function sample_balls(n::Int64, k::Int64, num_samples::Int64 )
  count = 0
  sample_size = 1
  #true_ratio = k/n
  sample_list = fill(0,num_samples)
  #println("i: ",0,"  true_rato: ",true_ratio)
  #abs_diff_sum = 0.0
  #estimate_sum = 0.0
  for i = 1:num_samples
    sample = rand(1:n, sample_size )
    count = length(findall( x->x<=k, sample ))
    sample_list[i] = count
    k -= count
    #if i % print_interval == 0
    #  println("i: ",i,"  count: ",count)
    #end
  end
  sample_list
end
  
# If hypothesis is true, result should be a list which approximates the list [p_i for i=0:(num_samples-1)]
function average_sample_lists( num_sample_lists::Int64, n::Int64, k::Int64, num_samples::Int64 )
  sum_sample_lists = fill(0,num_samples)
  for j = 1:num_sample_lists
    sl=sample_balls( n, k, num_samples )
    sum_sample_lists .+= sl
  end
  sum_sample_lists/num_sample_lists
end

# If hypothesis is true, result should be a list of floats which is close to zero.
function test_hypothesis( num_sample_lists::Int64, n::Int64, k::Int64, num_samples::Int64 )
  asl = average_sample_lists( num_sample_lists, n, k, num_samples )
  tpl = [ k/n*(1-1.0/n)^i for i = 0:(num_samples-1) ]  # Theoretical probability list
  (asl - tpl)./tpl
end
