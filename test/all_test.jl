using Test, SafeTestsets

begin
    @safetestset "triple_loop" begin include("triple_loop_test.jl") end
    @safetestset "geometry" begin include("geometry_test.jl") end
    @safetestset "binning" begin include("binning_test.jl") end
end