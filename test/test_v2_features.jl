# V2.0 features: étale diagonalization, multiplication matrix / trace / norm
# TEST blocks, transferGW, unstable degrees.
#
# Provenance: M2 TESTs 34, 36-39 of A1BrouwerDegrees v2.1, on the development
# branch of Macaulay2 (inline asserts) + M2 1.25.11 capture runs for literal
# references. Univariate unstable-degree Gram matrices compare LITERALLY
# (fixed monomial order, no basis call), unlike the multivariate stable degrees.
# The transfer tests take their inputs from M2 TESTs 40-41, but derive the
# expected values independently ("consistency"): trace forms by Newton's
# identities, the determinant identity det = det(T)^n N(det M), the signature at
# the real embeddings, the defining entries Tr(M[i, j] b_s b_t), and Hecke's
# restriction of scalars.
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

@testset "Transfer over QQ[x]/(x^5-x-1), M2 Test 40 input (consistency)" begin
    Rp, (x,) = polynomial_ring(QQ, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^5 - x - 1]))
    M = matrix(A, 3, 3,
        [A(1), A(3 * x^2 + 4 * x^4), A(8 * x^3 + 4),
         A(3 * x^2 + 4 * x^4), A(5), A(1),
         A(8 * x^3 + 4), A(1), A(7 * x^2 + 3 * x)])
    GQ = transfer_gw(GWClass(M))
    @test form_rank(GQ) == 15                     # rank 3 times [L:QQ] = 5
    # entry (5(i-1)+s, 5(j-1)+t) is Tr(M[i, j] b_s b_t) for b = 1, x, ..., x^4
    b = [A(x^m) for m in 0:4]
    G = gw_matrix(GQ)
    @test all(G[5 * (i - 1) + s, 5 * (j - 1) + t] == algebra_trace(A, M[i, j] * b[s] * b[t])
              for i in 1:3, j in 1:3, s in 1:5, t in 1:5)
    # trace form of L/QQ: Tr(x^m) = 5, 0, 0, 0, 4, 5, 0, 0, 4 for m = 0..8
    # (Newton's identities, then x^5 = x + 1)
    T = matrix(QQ, [5 0 0 0 4; 0 0 0 4 5; 0 0 4 5 0; 0 4 5 0 0; 4 5 0 0 4])
    @test det(G) == det(T)^3 * algebra_norm(A, det(M))
    # x^5 - x - 1 has one real root, at which M has signature 1; each pair of
    # complex embeddings contributes only hyperbolic planes
    @test form_signature(GQ) == 1
    # the transfer is the restriction of scalars of the quadratic space
    Rz, z = polynomial_ring(QQ, "z")
    K, r = number_field(z^5 - z - 1, "r")
    MK = map_entries(e -> evaluate(lift(e), [r]), M)
    V, _ = restrict_scalars(quadratic_space(K, MK), QQ)
    @test is_isomorphic_form(GQ, GWClass(gram_matrix(V)))
end

@testset "Transfer over GF(7)[x]/(x^3+6x^2+4), M2 Test 41 input (consistency)" begin
    k = GF(7)
    Rp, (x,) = polynomial_ring(k, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^3 + 6 * x^2 + 4]))   # irreducible: A is GF(7^3)
    M = matrix(A, 3, 3,
        [A(1), A(2), A(x),
         A(2), A(x^2 + 5), A(3 * x + 2),
         A(x), A(3 * x + 2), A(5)])
    G7 = transfer_gw(GWClass(M))
    @test form_rank(G7) == 9
    # trace form on 1, x, x^2: Tr(x^m) = 3, 1, 1, 3, 6 for m = 0..4 (Newton's
    # identities); its determinant 4 is a square, and N(det M) = 1
    T = matrix(k, [3 1 1; 1 1 3; 1 3 6])
    @test transfer_gw(GWClass(matrix(A, 1, 1, [A(1)]))) == GWClass(T)
    @test det(T) == k(4)
    @test algebra_norm(A, det(M)) == one(k)
    @test det(gw_matrix(G7)) == det(T)^3 * algebra_norm(A, det(M))
    # over a finite field, rank and discriminant classify forms
    @test is_isomorphic_form(G7, diagonal_form(k, (1, 1, 1, 1, 1, 1, 1, 1, 1)))
end

@testset "Transfer cross-checks (consistency)" begin
    Rp, (x,) = polynomial_ring(QQ, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^5 - x - 1]))
    # <1> goes to the trace form of A/QQ (Newton's identities, as above)
    T = matrix(QQ, [5 0 0 0 4; 0 0 0 4 5; 0 0 4 5 0; 0 4 5 0 0; 4 5 0 0 4])
    @test transfer_gw(GWClass(matrix(A, 1, 1, [A(1)]))) == GWClass(T)
    # the rank is multiplied by [A:QQ] = 5, and a diagonal form goes to the direct
    # sum of the transfers of its entries
    D = transfer_gw(GWClass(matrix(A, 2, 2, [A(x^4), A(0), A(0), A(x + 1)])))
    @test form_rank(D) == 10
    @test D == gw_direct_sum(transfer_gw(GWClass(matrix(A, 1, 1, [A(x^4)]))),
                             transfer_gw(GWClass(matrix(A, 1, 1, [A(x + 1)]))))
    # an entry of trace zero is no obstruction (Tr(x) = 0): a nondegenerate form
    # over a field transfers to a nondegenerate form
    Mx = matrix(A, 2, 2, [A(x), A(1), A(1), A(x^2)])
    @test form_rank(transfer_gw(GWClass(Mx))) == 10
    # well defined on classes: a congruent Gram matrix has an isomorphic transfer
    P = matrix(A, 2, 2, [A(1), A(x), A(0), A(1)])
    @test is_isomorphic_form(transfer_gw(GWClass(transpose(P) * Mx * P)),
                             transfer_gw(GWClass(Mx)))
    # <1> and <x^2> = <2> are the same class over QQ(sqrt 2), and so are their
    # transfers <2, 4> and <4, 8>
    B, _ = quo(Rp, ideal(Rp, [x^2 - 2]))
    t1 = transfer_gw(GWClass(matrix(B, 1, 1, [B(1)])))
    t2 = transfer_gw(GWClass(matrix(B, 1, 1, [B(x^2)])))
    @test t1 == diagonal_form(QQ, (2, 4))
    @test t2 == diagonal_form(QQ, (4, 8))
    @test is_isomorphic_form(t1, t2)
    # over QQ[t]/(t^2 - 1) = QQ x QQ, the sum of the specializations t = 1, t = -1
    St, (t,) = polynomial_ring(QQ, ["t"])
    E, _ = quo(St, ideal(St, [t^2 - 1]))
    beta = GWClass(matrix(E, 2, 2, [E(1), E(2), E(2), E(t)]))
    @test transfer_gw(beta) ==
        GWClass(matrix(QQ, [2 0 4 0; 0 2 0 4; 4 0 0 2; 0 4 2 0]))
    @test is_isomorphic_form(transfer_gw(beta),
        gw_direct_sum(GWClass(QQ[1 2; 2 1]), GWClass(QQ[1 2; 2 -1])))
    # 1 + t is a zero divisor (it maps to (2, 0)), so <1 + t> is degenerate on a
    # factor and has no transfer
    @test_throws ErrorException transfer_gw(GWClass(matrix(E, 1, 1, [E(1 + t)])))
end

@testset "Trace, norm and transfer over a quotient with lex ordering (consistency)" begin
    # QQ[x, y]/(x - y^2, y^4 - 2) is QQ(2^(1/4)); built with the lex ordering, its
    # standard monomials 1, y, y^2, y^3 are not those of the default ordering, and
    # the results must not depend on the ordering the quotient was built with
    S, (x, y) = polynomial_ring(QQ, ["x", "y"])
    I = ideal(S, [x - y^2, y^4 - 2])
    L, _ = quo(S, I; ordering = lex(S))
    @test algebra_trace(L, L(1)) == 4
    @test algebra_trace(L, L(x)) == 0              # Tr(2^(1/2))
    @test algebra_norm(L, L(y)) == -2              # N(2^(1/4)) = -2 from y^4 - 2
    Ld, _ = quo(S, I)
    @test multiplication_matrix(L, L(x + y)) == multiplication_matrix(Ld, Ld(x + y))
    # trace form on 1, x, y, xy, from Tr(y^m) = 4, 0, 0, 0, 8, 0, 0 for m = 0..6
    @test transfer_gw(GWClass(matrix(L, 1, 1, [one(L)]))) ==
        GWClass(matrix(QQ, [4 0 0 0; 0 8 0 0; 0 0 0 8; 0 0 8 0]))
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
    @test is_isomorphic_form(deg1, GWuClass(matrix(F, 1, 1, [inv(F(3))])))          # M2 TEST
    @test is_isomorphic_form(deg2, GWuClass(matrix(F, 1, 1, [F(8) * inv(F(3))])))   # M2 TEST
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
