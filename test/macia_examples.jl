# Use entroy.jl to reproduce the examples on page 397 of the Macia and Sole 2009 paper
# Note that this paper uses natural logarithms (logs to base exp(1).

using Test

const MyInt = UInt8

# Populations described by Table 2
X2_1 = [0x3,0x2,0x3,0x1]      # X^2_1
Xh2_1 = [0x1,0x0,0x1,0x1]     # \hat{X}^2_1
X = [0x7, 0x4, 0x7, 0x3]
O = [0x3, 0x3, 0x1, 0x0]

# Test the results given in Table 3 of Macia & Sole
@test entropy(X2_1,base=exp(1)) ≈ 1.0397207708399179
@test entropy(Xh2_1,base=exp(1)) ≈ 0.5623351446188083
@test entropy(X,base=exp(1)) ≈ 1.0397207708399179
@test entropy(O,base=exp(1)) ≈ 1.0397207708399179
@test joint_entropy(X,O,base=exp(1)) ≈ 1.3862943611198906 
@test mutual_information(X,O,base=exp(1)) ≈ 0.6931471805599452 
@test mutual_information(X2_1,O,base=exp(1)) + mutual_information(Xh2_1,O,base=exp(1)) - mutual_information(X,O,base=exp(1)) ≈ 0.2157615543388356 
