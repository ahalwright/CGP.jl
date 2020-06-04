#  Tests of the functions in Entropy.jl
using Test
using DataFrames
using Main.CGP

P = [1,1,1,2];
Q = [2,2,3,3];
PQ = [1,1,1,1,2,2,3,3];
if Main.CGP.MyInt == UInt8
  R = [0x01,0x01,0x01,0x01,0x02,0x02,0x03,0x03];
elseif Main.CGP.MyInt == UInt16
  R = [0x001,0x001,0x001,0x001,0x002,0x002,0x003,0x003];
end
S = ["A", "B", "C", "C"]
F1 = [0.2, 0.4, 0.1, 0.3]
F2 = [0.5, 0.2, 0.1, 0.2]
@test pop_to_dist(P) == Dict(1=>0.75,2=>0.25)  # Note 0.25 and 0.75 are represented exactly
@test pop_to_dist(Q) == Dict(2=>0.5,3=>0.5)
@test pop_to_dist(PQ) == Dict(1=>0.5,2=>0.25,3=>0.25)
@test pop_to_dist(S) == Dict("A"=>0.25, "B"=>0.25, "C"=> .5 )
@test pop_to_dist(R) == Dict( R[6]=>0.25, R[8]=>0.25, R[1]=> .5 )
@test pops_to_dist(P,Q) == Dict((1,2)=>0.5, (2,3)=>0.25, (1,3)=> .25 )
@test isapprox( pops_to_tbl([P,Q]), [0.375 0.125 0.0; 0.0 0.25 0.25] )
df = DataFrame(names1=["A","B","C"], counts1=[2,5,3], names2=["A","B","D"], counts2=[4,5,1])
@test isapprox( pop_counts_to_tbl(df), [0.1 0.25 0.15 0.0; 0.2 0.25 0.0 0.05])
tbl = [0.125 0.25 0.125 0.0; 0.125 0.25 0.0 0.125]
@test table_row_to_dist(tbl,1) == Dict(1=>0.25,2=>0.5,3=>0.25,4=>0.0)
@test column_marginal(tbl) == [0.25, 0.5, 0.125, 0.125]
@test isapprox( entropy(PQ), 1.5 )
tbl = [0.25 0.5 0.125 0.125; 0.375 0.25 0.0 0.125]
@test isapprox( entropy(tbl,1), 1.75 )
@test isapprox( entropy(tbl,2), 1.4591479170272448 )
@test relative_entropy(pop_to_dist(Q),pop_to_dist(PQ)) == 1.0
@test relative_entropy(Q,PQ) == 1.0
@test isapprox( entropy(F1), 1.8464393446710154 )
@test isapprox( entropy(F2), 1.7609640474436812 )
@test isapprox( mutual_information(F1,F2), 0.40563906222956625 )

# Test the Cover and Thomas (1991) implementations of joint, conditional, relative entropy; mutual information 
# Example 2.2.1 page 17 of Cover and Joy
tbl = [ 1/8 1/16 1/32 1/32; 1/16 1/8 1/32 1/32; 1/16 1/16 1/16 1/16; 1/4 0 0 0]
trans_tbl = convert(Array{Float64,2},transpose(tbl))
@test joint_entropy(tbl) ≈ 27/8  # H(X,Y)
X = column_marginal(tbl)
Y = row_marginal(tbl)
@test entropy(X) ≈ 7/4
@test entropy(Y) ≈ 2.0
@test conditional_entropy(tbl) ≈ 11/8   # H(X|Y)
@test conditional_entropy(trans_tbl) ≈ 13/8  # H(Y|X)
r = rand()
s = rand()
#p = Dict{Int64,Float64}( 0=>1.0-r, 1=>r )
#q = Dict{Int64,Float64}( 0=>1.0-q, 1=>q )
p = Main.CGP.DIST_TYPE( 0=>1.0-r, 1=>r )
q = Main.CGP.DIST_TYPE( 0=>1.0-s, 1=>s )
@test relative_entropy(p,q) ≈ (1-r)*log2((1-r)/(1-s)) + r*log2(r/s)  # Eq 2.31
@test relative_entropy(q,p) ≈ (1-s)*log2((1-s)/(1-r)) + s*log2(s/r)  # Eq 2.32
@test mutual_information(tbl) ≈ 3/8
@test mutual_information(trans_tbl) ≈ 3/8
@test mutual_information(tbl) ≈ entropy(X)-conditional_entropy(tbl)  # I(X;Y)=H(X)-H(X|Y) eq 2.43
@test mutual_information(tbl) ≈ entropy(Y)-conditional_entropy(trans_tbl)  # I(X;Y)=H(Y)-H(Y|X) eq 2.44
@test mutual_information(tbl) ≈ entropy(X)+entropy(Y)-joint_entropy(tbl)  # I(X;Y)=H(X)+H(Y)-H(X,Y0 eq 2.45
# Now test relationship of Cover mutual information to Sherwin mutual information
@test mutual_information(tbl) ≈ mutual_information(trans_tbl)  # mutual information is symmetric
@test mutual_information(tbl) ≈ sherwin_mutual_information(tbl)
@test mutual_information(trans_tbl) ≈ sherwin_mutual_information(trans_tbl)

