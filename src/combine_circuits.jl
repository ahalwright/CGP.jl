# Combines two 1-output circuits c1 and c2 with the same paramterers and number of gates into a 2-output circuit
# The output values of the result is the 2-output goal whose components are the output values of c1 and of c2.
function combine_circuits( c1::Chromosome, c2::Chromosome )
  # Two functions used only within this function
  inp1(i,numinputs) = i <= numinputs ? i : 2*(i-numinputs)+numinputs-1
  inp2(i,numinputs) = i <= numinputs ? i : 2*(i-numinputs)+numinputs
  @assert c1.params.numinputs == c2.params.numinputs
  @assert c1.params.numoutputs == 1
  @assert c2.params.numoutputs == 1
  if length(c1.interiors) < length(c2.interiors)
    println("  c1ints: ",length(c1.interiors),"  c2ints: ",length(c2.interiors))
    ctemp = c1
    c1 =c2
    c2 = ctemp
  end
  #println("c1ints: ",length(c1.interiors),"  c2ints: ",length(c2.interiors))
  p = Parameters( c1.params.numinputs, 2, c1.params.numinteriors+c2.params.numinteriors,
    c1.params.numlevelsback+c2.params.numlevelsback )
  inputs = c1.inputs
  interiors = InteriorNode[]
  new_numinteriors = 2*length(c2.interiors) 
  #println("new_numints: ",new_numinteriors)
  #print_parameters(p)
  for i = 1:length(c2.interiors)
    int_node1 = deepcopy(c1.interiors[i])
    input1 = inp1(c1.interiors[i].inputs[1],p.numinputs)
    input2 = inp1(c1.interiors[i].inputs[2],p.numinputs)
    int_node1.inputs = [input1,input2]
    #println("i: ",i,"  int_node1: ",int_node1)
    push!(interiors,int_node1)
    int_node2 = deepcopy(c2.interiors[i])
    input1 = inp2(c2.interiors[i].inputs[1],p.numinputs)
    input2 = inp2(c2.interiors[i].inputs[2],p.numinputs)
    int_node2.inputs = [input1,input2]
    #println("i: ",i,"  int_node2: ",int_node2)
    push!(interiors,int_node2)
  end
  #println("interiors: ",interiors)
  outputs = [OutputNode(p.numinputs+new_numinteriors-1),OutputNode(p.numinputs+new_numinteriors)]
  Chromosome( p, c1.inputs, interiors, outputs, 0.0, 0.0 )
end

