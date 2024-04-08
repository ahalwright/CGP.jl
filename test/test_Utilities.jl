using Test

# Import the function to be tested
include("/Users/oldmtnbiker/Library/CloudStorage/OneDrive-Personal/evotech/CGP.jl/src/Utilities.jl")

# Test the string_to_MyInt function
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