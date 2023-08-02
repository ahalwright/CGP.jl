# Attempts to compute a fitness non-decreasing path from the circuit represented by circ_int to a genotype of the phenotype phdest.
# Sometimes works, but often goes into a long loop.  
# However, our exact evolvability results from the evolvability paper show that there is a mutation between any pair of phenotypes
#    for the 3-input 8-gates case, which proves perfect navigability.
using DataStructures

function navigate( P::Parameters, funcs::Vector{Func}, circ_int::Int128, fitness::Vector{Float64}, phdest::Vector{MyInt} )
  @assert length(fitness) == 2^2^P.numinputs
  circ = int_to_chromosome( circ_int, P, funcs )
  phsrc = output_values( circ )
  fit = fitness[ phsrc[1]+1 ]
  println("phsrc: ",phsrc,"  src fit: ",fit)
  println("phdest: ",phdest,"  dest fit: ",fitness[ phdest[1]+1 ])
  path = Tuple{Int128,Float64}[(circ_int,fit)]
  queue = Queue{Vector{Tuple{Int128,Float64}}}()
  #explored = Set{Vector{Tuple{Int128,Float64}}}( [] )
  explored = Set{Int128}( [] )
  enqueue!( queue, path )
  while length(queue) > 0
    #println("queue: ", queue )
    current_path = dequeue!( queue )
    push!( explored, current_path[end][1] )
    current_circ = int_to_chromosome( current_path[end][1], P, funcs )
    #print("current_circ:  ")
    #print_circuit(current_circ )
    #println("prev fit: ",current_path[end][2])
    mut_circs = mutate_all( current_circ, funcs, output_outputs=false, output_circuits=true )
    for mcirc in mut_circs
      mcirc_int = chromosome_to_int( mcirc, funcs )
      if !(mcirc_int in explored)
        mph = output_values( mcirc )
        fit = fitness[ mph[1]+1 ]
        #print("mcirc:  fit: ",fit,"  mph: ",mph,"   ")
        #print_circuit(mcirc)
        if fitness[phsrc[1]+1] <= fit && fit <= fitness[phdest[1]+1] && fit >= current_path[end][2]
          mcirc_int = chromosome_to_int( mcirc, funcs )
          mpath = push!(path,(mcirc_int,fit))
          println("mph: ",mph,"  mcirc_int: ",mcirc_int,"  phdest: ",phdest)
          if mph == phdest
            return mpath
          end
          #if !(mpath in explored) 
          enqueue!( queue, mpath )
          push!( explored, mcirc_int )
          #end
        end
      end
    end
  end
end
