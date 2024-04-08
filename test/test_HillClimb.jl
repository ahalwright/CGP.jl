# test src/HillClimb.jl.   (also linked to src/hill_climbing.jl).
using Test

isapproxsigfigs(a, b, precision) = round(a, sigdigits=precision) == round(b, sigdigits=precision)
# Import the function to be tested
include("/Users/oldmtnbiker/Library/CloudStorage/OneDrive-Personal/evotech/CGP.jl/src/HillClimb.jl")
P2 = Parameters(2,1,6,3); funcs2 = default_funcs(P2)
P3 = Parameters(3,1,8,4); funcs3 = default_funcs(P3)
Random.seed!(1)
fitness = map( _->rand(), 1:2^2^P3.numinputs )

# Test the function hill_climb
@testset "HillClimb2" begin
Random.seed!(1)
rch2 = random_chromosome(P2,funcs2)
max_iterations = 20
(hch2, hfit2, hph2 ) = hill_climb( fitness, rch2, funcs2, max_iterations )
@test isapproxsigfigs(hfit2, 0.9149290036628314, 6)  
hcr2 = hill_climb_rank( fitness, rch2, funcs2, max_iterations )
@test hcr2 == 16
end

@testset "HillClimb3" begin
Random.seed!(1)
rch3 = random_chromosome(P3,funcs3)
max_iterations = 20
(hch3, hfit3, hph3 ) = hill_climb( fitness, rch3, funcs3, max_iterations )
@test isapproxsigfigs(hfit3, 0.9967010737682973, 6)  
hcr3 = hill_climb_rank( fitness, rch3, funcs3, max_iterations )
@test hcr3 == 256
end

#=  Example to be copied
@testset "string_to_MyInt tests" begin
    @test string_to_MyInt("10") == MyInt(10)
    @test string_to_MyInt("[20, 30, 40]") == MyInt(20)
    
    # Test for an illegal argument
    @test_throws ArgumentError string_to_MyInt("abc")
end

# Test the string_to_expression function
@testset "string_to_expression tests" begin
    @test string_to_expression("1 + 2") == 3
    @test string_to_expression("sqrt(4)") == 2.0
    
    # Test for an illegal argument
    @test_throws ArgumentError string_to_expression("abc")
end
=#
