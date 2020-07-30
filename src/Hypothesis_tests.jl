# Do two-sided equal variances T tests on data from a CSV file.
using CSV
using DataFrames
using HypothesisTests

# Returns a dafaframe by reading a csv file
# make sure that comment lines are preceded by #
function import_csv( file::String )
  CSV.read(file,comment="#")
end

function two_sided_t_test_pvalue( value_field::CSV.Column{Int64,Int64}, discriminator_field::CSV.Column{Bool,Bool} )
  pvalue( EqualVarianceTTest( value_field[ discriminator_field ], value_field[ .~discriminator_field ] ))
end

function two_sided_t_test( value_field::CSV.Column{Int64,Int64}, discriminator_field::CSV.Column{Bool,Bool} )
   EqualVarianceTTest( value_field[ discriminator_field ], value_field[ .~discriminator_field ] )
end    

function two_sided_t_test_pvalue( value_field::CSV.Column{Float64,Float64}, discriminator_field::CSV.Column{Bool,Bool} )
  pvalue( EqualVarianceTTest( value_field[ discriminator_field ], value_field[ .~discriminator_field ] ))
end

function two_sided_t_test( value_field::CSV.Column{Float64,Float64}, discriminator_field::CSV.Column{Bool,Bool} )
   EqualVarianceTTest( value_field[ discriminator_field ], value_field[ .~discriminator_field ] )
end    

# Do specific tests for nactive, complexity, degeneracy, and steps using fault_tol as the discriminator field
function nactive_test( file::String )
  df = import_csv( file )
  pvalue( EqualVarianceTTest( df.nactive[ df.fault_tol ], df.nactive[ .~df.fault_tol ] ) )
end

function complexity_test( file::String )
  df = import_csv( file )
  pvalue( EqualVarianceTTest( df.complexity[ df.fault_tol ], df.complexity[ .~df.fault_tol ] ) )
end

function degeneracy_test( file::String )
  df = import_csv( file )
  pvalue( EqualVarianceTTest( df.degeneracy[ df.fault_tol ], df.degeneracy[ .~df.fault_tol ] ) )
end

function steps_test( file::String )
  df = import_csv( file )
  pvalue( EqualVarianceTTest( df.steps[ df.fault_tol ], df.steps[ .~df.fault_tol ] ) )
end
