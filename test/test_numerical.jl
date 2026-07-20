# The ℂ numerical layer — only the UNSTABLE degrees use it. Over ℂ the stable
# degree is the identity form of the algebra dimension (a discrete invariant
# from the exact computation, tested in test_degrees.jl), so it needs no
# numerical path. The unstable degree carries a k^×-scalar and so must find the
# actual complex roots; those are the methods exercised here.
#
# Validation policy: NOT cross-validated against the exact core — computations
# over inexact fields are inexact. Tested against Macaulay2 1.25.11 CC runs up
# to small numerical tolerance (ranks exact, scalars |Δ| < 1e-8). Inputs are
# HomotopyContinuation.jl expressions.

using Test
using MotivicHomotopy
import HomotopyContinuation as HC
import LinearAlgebra

HC.@var x

@testset "CC unstable global degree (M2 capture, tolerance)" begin
    f = (x - 1) * (x - 2) * (x - 3)
    g = (x - 1) * (x - 4)
    G = global_unstable_A1_degree(f, g)
    @test gw_matrix(G) == Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)    # M2: rank 2
    @test abs(gw_scalar(G) - (-2.0)) < 1e-8       # M2: toCC(-.2p53e1, ~1e-64)

    fd = (x - 1)^2 * (x - 5)
    gd = x - 2
    GG = global_unstable_A1_degree(fd, gd)
    @test size(gw_matrix(GG), 1) == 3                                  # M2: rank 3
    @test abs(gw_scalar(GG) - (-2.9999999999999982)) < 1e-8            # M2 value
end

@testset "CC unstable local degrees (M2 capture, tolerance)" begin
    f = (x - 1) * (x - 2) * (x - 3)
    g = (x - 1) * (x - 4)
    d2 = local_unstable_A1_degree(f, g, 2)
    @test size(gw_matrix(d2), 1) == 1                                  # M2: rank 1
    @test abs(gw_scalar(d2) - 2.0) < 1e-8         # M2: toCC(.2p53e1, ~1e-64)

    fd = (x - 1)^2 * (x - 5)
    gd = x - 2
    L1 = local_unstable_A1_degree(fd, gd, 1)      # double root
    @test size(gw_matrix(L1), 1) == 2                                  # M2: rank 2
    @test abs(gw_scalar(L1) - (-0.0625)) < 1e-8   # M2: -1/16
    L5 = local_unstable_A1_degree(fd, gd, 5)
    @test size(gw_matrix(L5), 1) == 1                                  # M2: rank 1
    @test abs(gw_scalar(L5) - 0.1875) < 1e-8      # M2: 3/16
end

@testset "CC error paths (M2 UnstableLocalGlobalDegrees.m2 source)" begin
    f = (x - 1) * (x - 2) * (x - 3)
    g = (x - 1) * (x - 4)
    # not pointed after the common root cancels
    @test_throws ErrorException global_unstable_A1_degree(x - 1, (x - 1) * (x - 4))
    # r = 1 is a common root: cancelled, so not a zero after reduction
    @test_throws ErrorException local_unstable_A1_degree(f, g, 1)
    # linearTolerance must be positive
    @test_throws ErrorException local_unstable_A1_degree(f, g, 2; linear_tolerance = 0.0)
end
