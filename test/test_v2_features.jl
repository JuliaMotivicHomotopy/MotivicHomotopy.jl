# V2.0 features: étale diagonalization, multiplication matrix / trace / norm
# TEST blocks, transferGW, unstable degrees.
#
# Provenance: M2 TESTs 34, 36-41 (inline asserts) + M2 1.25.11 capture runs for
# literal references. Univariate unstable-degree Gram matrices compare LITERALLY
# (fixed monomial order, no basis call), unlike the multivariate stable degrees.
#
# Porting notes:
# - multiplication_matrix uses M2's basis order (ascending standard monomials), so
#   Tests 38/39 compare literally.
# - Test 38: M2's F = frac(GF(2)[x]/(x^2+x+1)) is the field GF(4); ported with
#   GF(2,2), whose generator satisfies the same minimal polynomial x^2+x+1.
# - Test 39: ported over the number field QQ[x]/(x^5-x-1); the expected values
#   depend only on the relation y^3 = -3y-2, not on x's minimal polynomial.

using Test
using Oscar
using MotivicHomotopy

@testset "M2 Test 34 — diagonalization over étale algebras" begin
    Rp, (x,) = polynomial_ring(QQ, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^2 + 1]))     # étale domain -> field path
    Sp, (y,) = polynomial_ring(QQ, ["y"])
    S, _ = quo(Sp, ideal(Sp, [y^2 - 1]))     # étale non-domain -> ring path
    G1 = GWClass(matrix(A, 2, 2, [A(1), A(2), A(2), A(x)]))
    G2 = GWClass(matrix(S, 2, 2, [S(1), S(2), S(2), S(y)]))
    G3 = GWuClass(matrix(A, 2, 2, [A(1), A(2), A(2), A(x)]))
    G4 = GWuClass(matrix(S, 2, 2, [S(1), S(2), S(2), S(y)]))
    DA = matrix(A, 2, 2, [A(1), A(0), A(0), A(x - 4)])
    DS = matrix(S, 2, 2, [S(1), S(0), S(0), S(y - 4)])
    @test diagonal_class(G1) == GWClass(DA)               # M2: === makeGWClass {{1,0},{0,x-4}}
    @test diagonal_class(G2) == GWClass(DS)
    @test diagonal_class(G3) == GWuClass(DA, gw_scalar(G3))
    @test diagonal_class(G4) == GWuClass(DS, gw_scalar(G4))
    @test stable_part(G3) == G1
    @test stable_part(G4) == G2
    # pivot != 1 distinguishes the two algorithms (M2 capture):
    # ring path (non-domain) is fraction-free; field path (domain) divides
    @test diagonalize_via_congruence(matrix(S, 2, 2, [S(2), S(1), S(1), S(y)])) ==
        matrix(S, 2, 2, [S(2), S(0), S(0), S(4 * y - 2)])
    @test diagonalize_via_congruence(matrix(A, 2, 2, [A(2), A(1), A(1), A(x)])) ==
        matrix(A, 2, 2, [A(2), A(0), A(0), A(x) - A(1) * inv(A(2))])
end

@testset "M2 Test 38 — multiplication matrix over the GF(4) tower" begin
    F4 = GF(2, 2)
    g = gen(F4)                                # g^2 + g + 1 == 0, M2's x
    @test iszero(g^2 + g + 1)
    Pk, (yv,) = polynomial_ring(F4, ["y"])
    K, _ = quo(Pk, ideal(Pk, [yv^2 + g * yv + 1]))
    N = multiplication_matrix(K, K(1 + g * yv))
    @test N == matrix(F4, 2, 2, [F4(1), g, g, 1 + g^2])   # M2: {{1, x}, {x, 1+x^2}}
end

@testset "M2 Test 39 — trace/norm over a quintic tower (domain substitute)" begin
    Ru, t = polynomial_ring(QQ, "t")
    NF, xnf = number_field(t^5 - t - 1, "x")   # see the Test 39 porting note above
    Py, (yv,) = polynomial_ring(NF, ["y"])
    K, _ = quo(Py, ideal(Py, [yv^3 + 3 * yv + 2]))
    a = K(1 + xnf * yv)
    N = multiplication_matrix(K, a)
    expected = matrix(NF, 3, 3,
        [NF(1), NF(0), -2 * xnf,
         xnf,   NF(1), -3 * xnf,
         NF(0), xnf,   NF(1)])                 # M2: {{1,0,-2x},{x,1,-3x},{0,x,1}}
    @test N == expected
    @test algebra_trace(K, a) == NF(3)         # M2: getTrace == 3_F
    @test algebra_norm(K, a) == det(N)         # M2: getNorm == det N
end

@testset "M2 Test 40 — transferGW over QQ[x]/(x^5-x-1)" begin
    Rp, (x,) = polynomial_ring(QQ, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^5 - x - 1]))
    M = matrix(A, 3, 3,
        [A(1), A(3 * x^2 + 4 * x^4), A(8 * x^3 + 4),
         A(3 * x^2 + 4 * x^4), A(5), A(1),
         A(8 * x^3 + 4), A(1), A(7 * x^2 + 3 * x)])
    GQ = transfer_gw(GWClass(M))
    # M2: === makeDiagonalForm(QQ, (5, -75, 17059280/279299)) — literal
    @test GQ == diagonal_form(QQ, (5, -75, QQ(17059280, 279299)))
end

@testset "M2 Test 41 — transferGW over GF(7)[x]/(x^3+6x^2+4)" begin
    k = GF(7)
    Rp, (x,) = polynomial_ring(k, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^3 + 6 * x^2 + 4]))
    M = matrix(A, 3, 3,
        [A(1), A(2), A(x),
         A(2), A(x^2 + 5), A(3 * x + 2),
         A(x), A(3 * x + 2), A(5)])
    G7 = transfer_gw(GWClass(M))
    @test is_isomorphic_form(G7, diagonal_form(k, (3, 4, 4)))
    @test gw_matrix(G7) == matrix(k, 3, 3,
        [k(3), k(0), k(0), k(0), k(-3), k(0), k(0), k(0), k(-3)])   # M2 capture
end

@testset "Transfer cross-checks (independent of M2)" begin
    Rp, (x,) = polynomial_ring(QQ, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^5 - x - 1]))
    # rank is preserved by M2's construction (entrywise trace of a diagonalization);
    # for x^5-x-1: tr(x^4) = 4, tr(x+1) = 5 (Newton's identities)
    M = matrix(A, 2, 2, [A(x^4), A(0), A(0), A(x + 1)])
    @test transfer_gw(GWClass(M)) == diagonal_form(QQ, (4, 5))
    @test form_rank(transfer_gw(GWClass(M))) == 2
    # trace of 1 in a degree-5 extension is 5: transfer of <1,1> is <5,5>
    I2 = matrix(A, 2, 2, [A(1), A(0), A(0), A(1)])
    @test transfer_gw(GWClass(I2)) == diagonal_form(QQ, (5, 5))
    # LIMITATION shared with M2: the entrywise trace can vanish on a nondegenerate
    # input, and the resulting degenerate diagonal form is rejected — M2's
    # makeGWClass errors identically (e.g. tr(x) = 0 here).
    Mbad = matrix(A, 2, 2, [A(x), A(1), A(1), A(x^2)])
    @test_throws ErrorException transfer_gw(GWClass(Mbad))
end

@testset "M2 Test 36 — unstable degrees over QQ (literal + class-level)" begin
    R, (x,) = polynomial_ring(QQ, ["x"])
    q = (x^5 - 6 * x^4 + 11 * x^3 - 2 * x^2 - 12 * x + 8) //
        (x^4 - 5 * x^2 + 7 * x + 1)
    Gdeg = global_unstable_A1_degree(q)
    M36 = matrix(QQ, [-68 38 11 -14 1; 38 -63 63 -29 7; 11 63 -84 39 -5;
                      -14 -29 39 -16 0; 1 7 -5 0 1])
    @test gw_matrix(Gdeg) == M36                    # M2 capture (literal)
    @test gw_scalar(Gdeg) == QQ(-53240)             # M2 capture
    @test is_isomorphic_form(Gdeg, GWuClass(M36, -53240))    # the M2 TEST assert
    deg1 = local_unstable_A1_degree(q, -1)
    deg2 = local_unstable_A1_degree(q, 1)
    deg3 = local_unstable_A1_degree(q, 2)
    @test gw_matrix(deg1) == matrix(QQ, [-5//27;;]) # M2 capture
    @test gw_matrix(deg2) == matrix(QQ, [-2;;])
    @test gw_matrix(deg3) == matrix(QQ, [0 0 11//3; 0 11//3 0; 11//3 0 0])
    @test gw_scalar(deg3) == QQ(-1331, 27)          # = det of the antidiagonal
    @test is_isomorphic_form(deg1, GWuClass(matrix(QQ, [-5//27;;])))   # M2 TEST
    @test is_isomorphic_form(deg2, GWuClass(matrix(QQ, [-2;;])))
    @test is_isomorphic_form(deg3, GWuClass(matrix(QQ, [0 0 11//3; 0 11//3 0; 11//3 0 0])))
    degSum = divisorial_sum([deg1, deg2, deg3], [-1, 1, 2])
    @test is_isomorphic_form(degSum, Gdeg)          # local-to-global (M2 TEST)
    @test gw_scalar(degSum) == QQ(-53240)           # M2 capture
end

@testset "M2 Test 37 — unstable degrees over GF(32003)" begin
    F = GF(32003)
    R, (x,) = polynomial_ring(F, ["x"])
    q = (x^2 + x - 2) // (3 * x + 5)
    Gdeg = global_unstable_A1_degree(q)
    @test gw_matrix(Gdeg) == matrix(F, [11 5; 5 3]) # M2 capture (literal)
    @test gw_scalar(Gdeg) == F(8)
    deg1 = local_unstable_A1_degree(q, -2)
    deg2 = local_unstable_A1_degree(q, 1)
    @test gw_matrix(deg1) == matrix(F, 1, 1, [inv(F(3))])        # M2: {{1/3}} = 10668
    @test gw_matrix(deg2) == matrix(F, 1, 1, [F(8) * inv(F(3))]) # M2: {{8/3}} = -10665
    @test is_isomorphic_form(Gdeg, GWuClass(matrix(F, [11 5; 5 3]))) # M2 TEST
    degSum = divisorial_sum([deg1, deg2], [-2, 1])
    @test is_isomorphic_form(degSum, Gdeg)
end

@testset "Stable part of unstable degrees = stable degrees (independent of M2)" begin
    R, (x,) = polynomial_ring(QQ, ["x"])
    f = x^3 - x
    @test is_isomorphic_form(stable_part(global_unstable_A1_degree(f)),
                             global_A1_degree([f]))
    # simple root
    @test is_isomorphic_form(stable_part(local_unstable_A1_degree(f, 1)),
                             local_A1_degree([f], ideal(R, [x - 1])))
    # multiple root: x^2 at the origin
    @test is_isomorphic_form(stable_part(local_unstable_A1_degree(x^2, 0)),
                             local_A1_degree([x^2], ideal(R, [x])))
    # same over GF(7)
    k = GF(7)
    S7, (s,) = polynomial_ring(k, ["s"])
    @test is_isomorphic_form(stable_part(global_unstable_A1_degree(s^2 - 1)),
                             global_A1_degree([s^2 - 1]))
end

@testset "Error paths (M2 UnstableLocalGlobalDegrees.m2 source)" begin
    R, (x,) = polynomial_ring(QQ, ["x"])
    # not pointed: deg f <= deg g
    @test_throws ErrorException global_unstable_A1_degree((x + 1) // (x^2 + 1))
    # local degree at a non-root
    @test_throws ErrorException local_unstable_A1_degree((x^2 - 1) // (x + 2), 3)
    # root not in the base field
    @test_throws ErrorException local_unstable_A1_degree(x^2 - 1, GF(5)(1))
end
