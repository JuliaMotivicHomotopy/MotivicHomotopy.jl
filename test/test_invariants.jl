# Invariants, in dependency order, plus the building-forms constructors that the
# M2 TEST blocks use to set up their inputs.
#
# Ported from the TEST blocks of the A1BrouwerDegrees Macaulay2 package (v2.0):
# Tests 7, 8, 16, 17, 21-25, 35, 42-46. Expected values are M2's asserts.
#
# Covered elsewhere: Tests 4-6 (test_degrees.jl); Tests 34, 36-37 and the étale
# paths (test_v2_features.jl); Test 47 + decomposition (test_decomposition.jl).

using Test
using Oscar
using MotivicHomotopy

# ---------- arithmetic ----------

@testset "M2 Test 16 — getPadicValuation" begin
    @test padic_valuation(27, 3) == 3
end

# ---------- rank / signature / discriminant / relevant primes ----------

@testset "M2 Test 44 — getRank" begin
    alpha = GWClass(matrix(QQ, [1 0 0; 0 -2 0; 0 0 5]))
    @test form_rank(alpha) == 3
    beta = GWClass([1.0 0; 0 -1.0])
    @test form_rank(beta) == 2
    gamma = GWClass(ComplexF64[1im 0 0 0; 0 -2 0 0; 0 0 3 0; 0 0 0 5])
    @test form_rank(gamma) == 4
    delta = GWClass(matrix(GF(7), [1 0; 0 3]))
    @test form_rank(delta) == 2
    epsilon = diagonal_form(QQ, (2, -3, 5, -7))    # M2: makeDiagonalForm
    @test form_rank(epsilon) == 4
    H = hyperbolic_form(QQ)                         # M2: makeHyperbolicForm
    @test form_rank(H) == 2
end

@testset "M2 Test 24 — signature, discriminant, relevant primes, Hasse-Witt" begin
    M1 = GWClass(matrix(QQ, [1 0 -3; 0 23 0; -3 0 -2//5]))
    M2 = GWClass(matrix(QQ, [1 0 0; 0 23 0; 0 0 -2//5]))
    M3 = GWClass(matrix(QQ, [1 0 0; 0 -23 0; 0 0 -2//5]))
    M4 = GWClass(matrix(QQ, [-1 0 0; 0 -23 0; 0 0 -2//5]))

    @test form_signature(M1) == 1
    @test form_signature(M2) == 1
    @test form_signature(M3) == -1
    @test form_signature(M4) == -3

    @test integral_discriminant(M1) == -5405
    @test relevant_primes(M1) == [23, 5, 47]
    @test hasse_witt_invariant(M1, 5) == -1
    @test hasse_witt_invariant(M1, 23) == 1
    @test hasse_witt_invariant(M1, 47) == -1
end

# ---------- Hilbert symbols ----------

@testset "M2 Test 25 — getHilbertSymbol / getHilbertSymbolReal" begin
    @test hilbert_symbol_padic(100, 7, 3) == 1
    @test hilbert_symbol_padic(100//121, 7//169, 3) == 1

    @test hilbert_symbol_padic(5, 1//9, 7) == 1
    @test hilbert_symbol_padic(1//9, 5, 7) == 1

    @test hilbert_symbol_padic(3, 11, 3) == -1
    @test hilbert_symbol_padic(3, 11, 2) == -1
    @test hilbert_symbol_padic(-3, -11, 2) == 1
    @test hilbert_symbol_padic(-5, 11, 2) == -1

    @test hilbert_symbol_real(-3//1, 5) == 1
    @test hilbert_symbol_real(-3, -5//1) == -1
    @test hilbert_symbol_real(-3//1, -5) == -1
    @test hilbert_symbol_real(3, 5) == 1
end

# ---------- anisotropic dimension ----------

@testset "M2 Test 45 — getAnisotropicDimension over QQ" begin
    H = hyperbolic_form(QQ)
    @test anisotropic_dimension(H) == 0
    twoPos = diagonal_form(QQ, (1, 1))
    @test anisotropic_dimension(twoPos) == 2
    hyp = diagonal_form(QQ, (1, -1))
    @test anisotropic_dimension(hyp) == 0
    fourPos = diagonal_form(QQ, (1, 1, 1, 1))
    @test anisotropic_dimension(fourPos) == 4
    M = GWClass(matrix(QQ, [1 0 0; 0 -2 0; 0 0 3]))
    @test anisotropic_dimension(M) + 2 * gw_witt_index(M) == form_rank(M)
end

@testset "M2 Test 46 — getAnisotropicDimensionQQp" begin
    twoPos = diagonal_form(QQ, (1, 1))
    @test anisotropic_dimension_qqp(twoPos, 2) == 2
    hyp = diagonal_form(QQ, (1, -1))
    @test anisotropic_dimension_qqp(hyp, 3) == 0
end

# ---------- Witt index (incl. the diagonalizeViaCongruence half of Test 17) ----------

@testset "M2 Test 17 — getWittIndex and diagonalizeViaCongruence" begin
    B = matrix(QQ, [0 1; 1 0])
    beta = GWClass(B)
    @test gw_witt_index(beta) == 1
    P = matrix(QQ, [0 5 1; 2 2 1; 0 0 1])
    A = matrix(QQ, [1 0 0; 0 -1 0; 0 0 1])
    @test gw_witt_index(GWClass(diagonalize_via_congruence(P * A * transpose(P)))) == 1
end

# ---------- isotropy ----------

@testset "M2 Test 21 — isIsotropic / isAnisotropic" begin
    A1 = matrix(QQ, [0 1; 1 0])
    @test is_isotropic_form(A1)
    @test !is_anisotropic_form(GWClass(A1))

    A2 = [1.0 -2 4; -2 2 0; 4 0 -7]
    @test !is_anisotropic_form(A2)
    @test is_isotropic_form(GWClass(A2))

    k = GF(13, 4)                                   # M2: GF(13^4)
    A3 = matrix(k, [7 81 63; 81 7 55; 63 55 109])
    @test is_isotropic_form(GWClass(A3))            # Chevalley-Warning

    A4 = matrix(QQ, [5 0; 0 5])
    @test is_anisotropic_form(A4)

    A5 = ComplexF64[3+1im 0; 0 5-1im]
    @test !is_anisotropic_form(A5)
end

# ---------- isomorphism (capstone) ----------

@testset "M2 Test 22 — isIsomorphicForm over QQ" begin
    B1 = matrix(QQ, [1 -2 4; -2 2 0; 4 0 -7])
    B2 = matrix(QQ, [-17198//4225 -166126//975 -71771//1560;
                     -166126//975 -27758641//4050 -251077//135;
                     -71771//1560 -251077//135 -290407//576])
    @test is_isomorphic_form(GWClass(B1), GWClass(B2))
    B3 = matrix(QQ, [-38 -50 23; -50 -62 41; 23 41 29])
    @test is_isomorphic_form(GWClass(B1), GWClass(B3))
end

@testset "M2 Test 23 — isIsomorphicForm over RR, QQ, GF(13)" begin
    A1 = [1.0 -2 4; -2 2 0; 4 0 -7]
    A2 = [-38.0 -50 23; -50 -62 41; 23 41 29]
    @test is_isomorphic_form(GWClass(A1), GWClass(A2))

    B1 = matrix(QQ, [1 -2 4; -2 2 0; 4 0 -7])
    B2 = matrix(QQ, [-38 -50 23; -50 -62 41; 23 41 29])
    @test is_isomorphic_form(GWClass(B1), GWClass(B2))

    k = GF(13)
    C1 = matrix(k, [1 11 4; 11 2 0; 4 0 6])
    C2 = matrix(k, [1 2 10; 2 3 2; 10 2 3])
    @test is_isomorphic_form(GWClass(C1), GWClass(C2))
end

@testset "M2 Test 7 — Pfister form vs 2H over GF(17)" begin
    P = pfister_form(GF(17), (2, 3))
    twoH = hyperbolic_form(GF(17), 4)
    @test is_isomorphic_form(P, twoH)
end

@testset "M2 Test 8 — diagonal and hyperbolic forms over RR" begin
    alpha = diagonal_form(Float64, (1, -1))   # M2: makeDiagonalForm(RR, (1,-1))
    beta = GWClass([0.0 1; 1 0])
    H = hyperbolic_form(Float64)              # M2: makeHyperbolicForm RR
    @test is_isomorphic_form(alpha, H)
    @test is_isomorphic_form(beta, H)
end

@testset "M2 Test 35 — isIsomorphicForm for unstable classes over QQ" begin
    B1 = matrix(QQ, [1 -2 4; -2 2 0; 4 0 -7])
    B2 = matrix(QQ, [-17198//4225 -166126//975 -71771//1560;
                     -166126//975 -27758641//4050 -251077//135;
                     -71771//1560 -251077//135 -290407//576])
    @test is_isomorphic_form(GWuClass(B1), GWuClass(B2, -18))
    B3 = matrix(QQ, [-38 -50 23; -50 -62 41; 23 41 29])
    @test is_isomorphic_form(GWuClass(B1), GWuClass(B3, -18))
end

@testset "M2 Test 42 — unstable builders + isIsomorphicForm over RR" begin
    alpha = diagonal_unstable_form(Float64, (1, -1))
    beta = GWuClass([0.0 1; 1 0])
    H = hyperbolic_unstable_form(Float64)
    @test is_isomorphic_form(alpha, H)
    @test is_isomorphic_form(beta, H)
end

@testset "M2 Test 43 — unstable builders + isIsomorphicForm over GF(27)" begin
    k = GF(3, 3)                               # M2: GF(27)
    alpha = diagonal_unstable_form(k, (1, -1))
    beta = GWuClass(matrix(k, [0 1; 1 0]))
    H = hyperbolic_unstable_form(k)
    @test is_isomorphic_form(alpha, H)
    @test is_isomorphic_form(beta, H)
end

@testset "Error paths (M2 invariants source)" begin
    # getIntegralDiscriminant / getRelevantPrimes / getSignature restrictions
    gf = GWClass(matrix(GF(7), [1 0; 0 3]))
    @test_throws ErrorException integral_discriminant(gf)
    @test_throws ErrorException relevant_primes(gf)
    @test_throws ErrorException form_signature(gf)
    # getHilbertSymbol validation
    @test_throws ErrorException hilbert_symbol_padic(0, 3, 5)
    @test_throws ErrorException hilbert_symbol_padic(2, 3, 6)
    @test_throws ErrorException hilbert_symbol_real(0, 3)
    # makeHyperbolicForm rejects odd rank
    @test_throws ErrorException hyperbolic_form(QQ, 3)
    # isIsomorphicForm: mismatched base fields
    @test_throws ErrorException is_isomorphic_form(GWClass(matrix(QQ, [1 0; 0 1])),
                                                   GWClass([1.0 0; 0 1]))
    # GWu linearTolerance must be positive (M2 validates in the unstable variant)
    a = GWuClass(matrix(QQ, [1 0; 0 1]))
    @test_throws ErrorException is_isomorphic_form(a, a; linear_tolerance = 0.0)
end
