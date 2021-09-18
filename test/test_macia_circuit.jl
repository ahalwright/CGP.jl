# Example of the circuit example of Figure 2 of the Macia & Sole (2009) paper.
#include("../../CGP.jl/src/CGP.jl")
using Main.CGP
#include("../../information_theory/src/entropy.jl")
using Combinatorics
using Test

# According to Macia & Sole's definition of Z as the number of interacting units,
#  you would think that Z = 5.  However, their Figure 2 example clearly states 
#  that X consists of gates 1, 2, and 3 (which are int3, int4, and int5 below).
#  Thus, Z = 3 even though gates 4 and 5 are computation nodes.
# As of March 13, I am using n from Tononi rather than Z from Macia & Sole
# This means results for degeneracy and complexity may not match Macia and Sole

# Runs without errors on 4/17/21
function test_macia_circuit()
  p0 = Parameters( numinputs=2, numoutputs=2, numinteriors=5, numlevelsback=8 )
  default_funcs( p0.numinputs )  # need to call default_funcs() to define global vars ONE and Ones
  # Circuit of Figure 2 of the Macia & Sole paper
  println("Figure 2 circuit")
  Ones = Main.CGP.construct_ones(p0.numinputs)[p0.numinputs]
  One() = Ones
  funcs = [ Func(Main.CGP.Nand, 2, "NAND") ]
  in1 = InputNode(1)
  in2 = InputNode(2)
  int3 = InteriorNode(NAND, [1, 2])
  int4 = InteriorNode(NAND, [3, 2])
  int5 = InteriorNode(NAND, [3, 2])
  int6 = InteriorNode(NAND, [4, 1])
  int7 = InteriorNode(NAND, [2, 5])
  out8 = OutputNode(6)
  out9 = OutputNode(7)
  c0 = Chromosome(p0,[in1,in2],[int3,int4,int5,int6,int7],[out8,out9],0.0,0.0)
  context = construct_contexts(p0.numinputs)[p0.numinputs]
  execute_chromosome(c0,context)
  println("Macia Sole Figure 2 chromosome: ")
  print_chromosome(c0)
  println("num_mutate_locations: ",num_mutate_locations(c0,funcs))
  @assert number_active(c0) == number_active_old(c0)
  nv = node_values(c0)
  println("nv: ",nv)
  
  # Node labels are from Macia Figure 2
  I1 = context[1]
  I2 = context[2]
  G1 = Nand(I1,I2)
  G2 = Nand(G1,I2)
  G3 = Nand(G1,I2)
  G4 = Nand(G2,I1)
  G5 = Nand(I2,G3)
  @test nv == ([I1,I2],[G1,G2,G3,G4,G5],[G4,G5])
  IN = nv[1];  
  XX = nv[2]   # 5 gates, while Macia uses only the first 3 gates for X
  X = XX[1:3]  # 3 gates, what Macia calls X
  O = nv[3]
  @test nv == (IN,XX,O)  # Failed on 8/15/21 because node caching in evaluate_node() is commented out
  println("Verifying the results in comparison to Tables 3 and 4 of Macia and Sole")
  gbX = get_bits( X, p0.numinputs )
  gbO = get_bits( O, p0.numinputs )
  X2_1 = X[1:2]
  Xh2_1 = X[3:3]
  gbX2_1 = get_bits( X2_1, p0.numinputs )
  gbXh2_1 = get_bits( Xh2_1, p0.numinputs )
  println("(gbX,gbO): ",(gbX,gbO))
  miX_O = mutual_information(gbX,gbO,base=exp(1))
  jeX_O = joint_entropy(gbX,gbO,base=exp(1))
  print("je(X,O): ",jeX_O)
  println("   mi(X,O): ",miX_O)
  println("(gbX2_1,gbO): ",(gbX2_1,gbO))
  jeX2_1_O = joint_entropy(gbX2_1,gbO,base=exp(1))
  miX2_1_O = mutual_information(gbX2_1,gbO,base=exp(1))
  println("je(X2_1,O): ",jeX2_1_O,"  mi(X2_1,O): ",miX2_1_O)
  jeXh2_1_O = joint_entropy(gbXh2_1,gbO,base=exp(1))
  miXh2_1_O = mutual_information(gbXh2_1,gbO,base=exp(1))
  println("je(Xh2_1,O): ",jeXh2_1_O,"  mi(Xh2_1,O): ",miXh2_1_O)
  println("je(X2_1,O) + je(Xh2_1,O) - je(X,O): ",jeX2_1_O + jeXh2_1_O - jeX_O)
  println("mi(X2_1,O) + mi(Xh2_1,O) - mi(X,O): ",miX2_1_O + miXh2_1_O - miX_O)
  #mi0 = map(x->(x,mutual_information(x,gbO)),map(x->get_bits(x,p0.numinputs),[s for s in combinations(X,2)]))
  #println("mi0: ",mi0)
  
  # Simplified and modified version of the above circuit with only one output gate
  println("Alternate simple chromosome: ")
  p1 = Parameters(1, 4, 0.05, 0.0, 2, 1, 2, 4, 8)
  in1 = InputNode(1)
  in2 = InputNode(2)
  int3 = InteriorNode(NAND, [1, 2])  # (in1,in2)
  int4 = InteriorNode(NAND, [3, 2])  # (int3,in2)
  int5 = InteriorNode(NAND, [1, 4])  # (in1,int4)
  int6 = InteriorNode(NAND, [5, 3])  # (int5,int3)
  out7 = OutputNode(6)
  c1 = Chromosome(p1,[in1,in2],[int3,int4,int5,int6],[out7],0.0,0.0)
  context = construct_contexts(p1.numinputs)[p1.numinputs]
  execute_chromosome(c1,context)
  #println("c1: ",c1)
  #print_chromosome(c1)
  println("num_mutate_locations: ",num_mutate_locations(c1,funcs))
  @assert number_active(c1) == number_active_old(c1)
  nv = node_values(c1)
  println("nv: ",nv)
  
  I1 = context[1]
  I2 = context[2]
  INT3 = Nand(I1,I2)
  INT4 = Nand(INT3,I2)
  INT5 = Nand(I1,INT4)
  INT6 = Nand(INT5,INT3)
  OUT7 = INT6
  @test nv == ([I1,I2],[INT3,INT4,INT5,INT6],[OUT7])
  IN = nv[1];  
  X = nv[2]
  O = nv[3]
  @test nv == (IN,X,O)
  gbX = get_bits( X, p1.numinputs )
  gbO = get_bits( O, p1.numinputs )
  println("(gbX,gbO): ",(gbX,gbO))
  println("mi(X,O): ",mutual_information(gbX,gbO))
  mi1 = map(x->(x,mutual_information(x,gbO)),map(x->get_bits(x,p1.numinputs),[s for s in combinations(X,2)]))
  println("mi1: ",mi1)
  
  # Duplicate the above circuit with each part going to a separate output
  p2 = Parameters(1, 4, 0.05, 0.0, 2, 2, 2, 8, 16)
  print_parameters(p2)
  in1 = InputNode(1)
  in2 = InputNode(2)
  int3 = InteriorNode(NAND, [1, 2])  # (in1,in2)
  int4 = InteriorNode(NAND, [3, 2])  # (int3,in2)
  int5 = InteriorNode(NAND, [1, 4])  # (in1,int4)
  int6 = InteriorNode(NAND, [5, 3])  # (int5,int3)
  int7 = InteriorNode(NAND, [1, 2])  # (in1,in2)
  int8 = InteriorNode(NAND, [7, 2])  # (int7,in2)
  int9 = InteriorNode(NAND, [1, 8])  # (in1,int8)
  int10 = InteriorNode(NAND, [9, 7])  # (int9,int7)
  out11 = OutputNode(6)
  out12 = OutputNode(10)
  c2 = Chromosome(p2,[in1,in2],[int3,int4,int5,int6,int7,int8,int9,int10],[out11,out12],0.0,0.0)
  context = construct_contexts(p2.numinputs)[p2.numinputs]
  execute_chromosome(c2,context)
  #println("c2: ",c2)
  print_chromosome(c2)
  println("num_mutate_locations: ",num_mutate_locations(c2,funcs))
  @assert number_active(c2) == number_active_old(c2)
  nv = node_values(c2)
  println("nv: ",nv)
  
  I1 = context[1]
  I2 = context[2]
  INT3 = Nand(I1,I2)
  INT4 = Nand(INT3,I2)
  INT5 = Nand(I1,INT4)
  INT6 = Nand(INT5,INT3)
  INT7 = Nand(I1,I2)
  INT8 = Nand(INT7,I2)
  INT9 = Nand(I1,INT8)
  INT10 = Nand(INT9,INT7)
  OUT11 = INT6
  OUT12 = INT10
  @test nv == ([I1,I2],[INT3,INT4,INT5,INT6,INT7,INT8,INT9,INT10],[OUT11,OUT12])
  IN = nv[1];  
  X = nv[2]
  O = nv[3]
  @test nv == (IN,X,O)
  gbX = get_bits( X, p2.numinputs )
  gbO = get_bits( O, p2.numinputs )
  gbX2_1 = get_bits(X2_1, p2.numinputs )
  println("(gbX,gbO): ",(gbX,gbO))
  println("mi(X,O): ",mutual_information(gbX,gbO))
  mi2 = map(x->(x,mutual_information(x,gbO)),map(x->get_bits(x,p2.numinputs),[s for s in combinations(X,2)]))
  println("mi2: ",mi2)
end
