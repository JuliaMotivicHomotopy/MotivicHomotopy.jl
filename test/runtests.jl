using Test

@testset "MotivicHomotopy.jl" begin
    include("test_core_types.jl")
    include("test_matrix_layer.jl")
    include("test_invariants.jl")
    include("test_decomposition.jl")
    include("test_degrees.jl")
    include("test_v2_features.jl")
    include("test_numerical.jl")
end
