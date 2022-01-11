# Test the 3-input 1-output case of Hu (2020)
# To run:  from evotech/CGP.jl/src:
# julia -L CGP.jl -L Fnc.jl ../test/testFnc.jl
using Random
using Test
using Main.CGP
@testset "Testing Fnc.jl for use_lincircuit==true and use_lincircuit==false, and both dict_to_csv() and dict_csv()" begin

p = Parameters(3,1,3,2) ;funcs=default_funcs(p)
pheno = 0x0023
Random.seed!(1); 
dtc_df = component_properties( p, pheno, funcs, use_lincircuit=true, nwalks_per_circuit=3, nwalks_per_set=3, walk_length=50, use_dict_csv=false )
println(dtc_df)
Random.seed!(1); 
dc_df = component_properties( p, pheno, funcs, use_lincircuit=true, nwalks_per_circuit=3, nwalks_per_set=3, walk_length=50, use_dict_csv=true )
println(dc_df)
#@test dc_df == dtc_df
pdf = pheno_counts_lc( p, funcs )
@test dataframe_count_phenos(dtc_df) == Float64(pdf[pdf.goal.==@sprintf("0x%04x",pheno),:counts][1])

p = Parameters(3,1,4,3) ;funcs=default_funcs(p)
pheno = pheno
Random.seed!(1); 
dtc_df = component_properties( p, pheno, funcs, use_lincircuit=false, nwalks_per_circuit=3, nwalks_per_set=3, walk_length=50, use_dict_csv=false )
println(dtc_df)
Random.seed!(1); 
dc_df = component_properties( p, pheno, funcs, use_lincircuit=false, nwalks_per_circuit=3, nwalks_per_set=3, walk_length=50, use_dict_csv=false )
println(dc_df)
#@test dc_df == dtc_df
pdf = pheno_counts_ch( p, funcs )
@test dataframe_count_phenos(dtc_df) == Float64(pdf[pdf.goal.==@sprintf("0x%04x",pheno),:counts][1])

end #testset
