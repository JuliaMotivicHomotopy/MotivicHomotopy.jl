# Matrix layer: diagonalization via congruence, per-field simplification,
# diagonal class/entries.
#
# Ported from the TEST blocks of the A1BrouwerDegrees Macaulay2 package (v2.0):
# Tests 0-2, 10-15, and the getDiagonalClass halves of Tests 26-28 (whose
# stable-part halves are in test_core_types.jl). Expected values are M2's asserts.
#
# Covered elsewhere:
# - Test 17: its diagonalizeViaCongruence assertion is wrapped in getWittIndex
#   (test_invariants.jl).
# - Test 34 and the diagonalizeViaCongruence doc example: étale algebras
#   (test_v2_features.jl).

using Test
using Oscar
using MotivicHomotopy
import LinearAlgebra

@testset "M2 Test 0 — getDiagonalClass over RR" begin
    M1 = [0.0 1; 1 0]
    G2 = diagonal_class(GWClass(M1))
    @test gw_matrix(G2) == [1.0 0; 0 -1]
end

@testset "M2 Test 1 — getDiagonalClass over CC" begin
    M3 = ComplexF64[1 2 3; 2 4 5; 3 5 7]
    G4 = diagonal_class(GWClass(M3))
    @test gw_matrix(G4) == ComplexF64[1 0 0; 0 1 0; 0 0 1]
end

@testset "M2 Test 2 — getDiagonalClass over QQ" begin
    M5 = matrix(QQ, [1 2 3; 2 4 5; 3 5 7])
    G6 = diagonal_class(GWClass(M5))
    @test gw_matrix(G6) == matrix(QQ, [1 0 0; 0 -2 0; 0 0 2])
end

@testset "M2 Test 10 — diagonal class/entries over CC" begin
    M1 = ComplexF64[1 0 0; 0 2 0; 0 0 -3]
    G = GWClass(M1)
    @test gw_matrix(diagonal_class(G)) == ComplexF64[1 0 0; 0 1 0; 0 0 1]
    @test diagonal_entries(G) == [1, 2, -3]
end

@testset "M2 Test 11 — diagonal class/entries over RR" begin
    M1 = [1.0 0 0; 0 2 0; 0 0 -3]
    G = GWClass(M1)
    @test gw_matrix(diagonal_class(G)) == [1.0 0 0; 0 1 0; 0 0 -1]
    @test diagonal_entries(G) == [1, 2, -3]
end

@testset "M2 Test 12 — diagonal class/entries over QQ" begin
    M = matrix(QQ, [1 0 0; 0 2 0; 0 0 -3])
    G = GWClass(M)
    @test gw_matrix(diagonal_class(G)) == M
    @test diagonal_entries(G) == [1, 2, -3]
end

@testset "M2 Test 13 — diagonal class/entries over GF(5)" begin
    # In GF(5), -1 is a square, so the nonsquare representative is taken from the
    # diagonal itself; diag(1,2,-3) is fixed by the simplification.
    M = matrix(GF(5), [1 0 0; 0 2 0; 0 0 -3])
    G = GWClass(M)
    @test gw_matrix(diagonal_class(G)) == M
    @test diagonal_entries(G) == [1, 2, -3]
end

@testset "M2 Test 14 — diagonal class/entries over GF(7)" begin
    kk = GF(7)
    M1 = matrix(kk, [1 0 0; 0 2 0; 0 0 -3])
    M2 = matrix(kk, [1 0 0; 0 1 0; 0 0 1])
    G = GWClass(M1)
    @test gw_matrix(diagonal_class(G)) == M2   # 1, 2, -3=4 are all squares mod 7
    @test diagonal_entries(G) == [1, 2, -3]
end

@testset "M2 Test 15 — squarefree reduction over QQ" begin
    M1 = matrix(QQ, [18 0 0; 0 125//9 0; 0 0 -8//75])
    M2 = matrix(QQ, [2 0 0; 0 5 0; 0 0 -6])
    G1 = GWClass(M1)
    @test gw_matrix(diagonal_class(G1)) == M2
end

@testset "M2 Test 26 — GWu getDiagonalClass over RR" begin
    M1 = [0.0 1; 1 0]
    G1 = GWuClass(M1)
    G2 = diagonal_class(G1)
    @test gw_matrix(G2) == [1.0 0; 0 -1]
    @test gw_scalar(G2) == LinearAlgebra.det(M1)   # M2: getScalar(G2) === det M1
end

@testset "M2 Test 27 — GWu getDiagonalClass over CC" begin
    M3 = ComplexF64[1 2 3; 2 4 5; 3 5 7]
    G3 = GWuClass(M3)
    G4 = diagonal_class(G3)
    @test gw_matrix(G4) == ComplexF64[1 0 0; 0 1 0; 0 0 1]
    @test gw_scalar(G4) == LinearAlgebra.det(M3)
end

@testset "M2 Test 28 — GWu getDiagonalClass over QQ" begin
    M5 = matrix(QQ, [1 2 3; 2 4 5; 3 5 7])
    G5 = GWuClass(M5)
    G6 = diagonal_class(G5)
    @test gw_matrix(G6) == matrix(QQ, [1 0 0; 0 -2 0; 0 0 2])
    @test gw_scalar(G6) == det(M5)
end

@testset "Diagonal-class cache (M2 documented behavior)" begin
    # M2 caches the diagonal class on beta.cache (GrothendieckWittClassesDoc.m2
    # shows beta.cache.getDiagonalClass). One deliberate difference: M2's first
    # call returns a fresh twin of the cached object; ours returns the cached
    # object itself (indistinguishable under ==).
    beta = GWClass(matrix(QQ, [1 2 3; 2 4 5; 3 5 7]))
    d1 = diagonal_class(beta)
    @test haskey(beta.cache, :diagonal_class)
    @test diagonal_class(beta) === d1
    betau = GWuClass(matrix(QQ, [1 2 3; 2 4 5; 3 5 7]))
    du = diagonal_class(betau)
    @test haskey(betau.cache, :diagonal_class)
    @test diagonal_class(betau) === du
end

@testset "Error paths (M2 diagonalizeViaCongruence source)" begin
    # M2: "matrix is not symmetric" for non-symmetric input
    @test_throws ErrorException diagonalize_via_congruence(matrix(QQ, [1 2; 3 4]))
    # (Étale diagonalization is covered by the M2 Test 34 port, test_v2_features.jl.)
end
