# Use Entropy.jl to reproduce the examples on page 397 of the Macia and Sole 2009 paper
# Note that this paper uses natural logarithms (logs to base exp(1).

using Test

p0 = Parameters( numinputs=2, numoutputs=2, numinteriors=5, numlevelsback=8 ) 
c0 = build_chromosome(    # Macia Sole Figure 2 chromosome
  (1,2),  ((NAND,Integer[1, 2]),(NAND,Integer[3, 2]),(NAND,Integer[3, 2]),(NAND,Integer[4, 1]),(NAND,Integer[2, 5])),  (6,7), 0.0)
out0 = output_values(c0)   # Executes c0
(IN,XX,O) = node_values(c0)
X = XX[1:3]  # the 3 gates that Macia calls X
X2_1 = X[1:2]
Xh2_1 = X[3:3]
gbX = get_bits( X, p0.numinputs )
gbO = get_bits( O, p0.numinputs )
gbX2_1 = get_bits( X2_1, p0.numinputs )
gbXh2_1 = get_bits( Xh2_1, p0.numinputs )          
@test gbX == MyInt[0x7, 0x4, 0x7, 0x3]      # gbX from ../test/macia_cicuit.jl
@test gbO == MyInt[0x3, 0x3, 0x1, 0x0] 
@test gbX2_1 == MyInt[0x3,0x2,0x3,0x1] 
@test gbXh2_1 == MyInt[0x1,0x0,0x1,0x1]
# Populations described by Table 2 

# Test the results given in Tables 3 and 4 of Macia & Sole
@test entropy(gbX2_1,base=exp(1)) ≈ 1.0397207708399179
@test entropy(gbXh2_1,base=exp(1)) ≈ 0.5623351446188083
@test entropy(gbX,base=exp(1)) ≈ 1.0397207708399179
@test entropy(gbO,base=exp(1)) ≈ 1.0397207708399179
@test joint_entropy(gbX,gbO,base=exp(1)) ≈ 1.3862943611198906 
@test mutual_information(gbX,gbO,base=exp(1)) ≈ 0.6931471805599452 
@test mutual_information(gbX2_1,gbO,base=exp(1)) + mutual_information(gbXh2_1,gbO,base=exp(1)) - mutual_information(gbX,gbO,base=exp(1)) ≈ 0.2157615543388356 
