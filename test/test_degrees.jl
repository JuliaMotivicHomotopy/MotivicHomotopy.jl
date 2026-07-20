# Local and global A¹-degrees via the Bézoutian (exact path: QQ, F_q).
#
# Ported from M2 TEST blocks 4, 5, 6, 9; Gram-matrix reference values captured from
# M2 1.25.11. Comparisons against captures are class-level (is_isomorphic_form):
# the Gram matrix depends on the standard-monomial basis ORDER, and Oscar's
# monomial_basis orders differently than M2's basis — the results are
# permutation-congruent. M2's own TEST blocks compare class-level too.
#
# ℝ/ℂ degree paths are excluded per the field model (ℝ: compute over QQ and
# base-change; ℂ: numerical layer, test_numerical.jl).

using Test
using Oscar
using MotivicHomotopy

@testset "M2 Test 4 — global degree of x^2 over QQ" begin
    T1, (x,) = polynomial_ring(QQ, ["x"])
    alpha = global_A1_degree([x^2])
    beta = GWClass(matrix(QQ, [0 1; 1 0]))    # M2 capture: exactly {{0,1},{1,0}}
    @test form_rank(alpha) == 2
    @test is_isomorphic_form(alpha, beta)
end

@testset "M2 Test 5 — two-variable system over QQ, local-global" begin
    T2, (z1, z2) = polynomial_ring(QQ, ["z1", "z2"])
    f1 = [(z1 - 1) * z1 * z2, QQ(3, 5) * z1^2 - QQ(17, 3) * z2^2]
    f1GD = global_A1_degree(f1)
    @test gw_witt_index(f1GD) == 3
    @test form_rank(f1GD) == 6                 # M2 capture
    I1 = ideal(T2, [z1, z2])
    I2 = ideal(T2, [z1 - 1, z2^2 - QQ(9, 85)])
    f1LD1 = local_A1_degree(f1, I1)
    f1LD2 = local_A1_degree(f1, I2)
    @test is_isomorphic_form(gw_direct_sum(f1LD1, f1LD2), f1GD)
    # M2 captures (class-level):
    @test is_isomorphic_form(f1LD1, GWClass(matrix(QQ,
        [0 0 0 17//3; 0 3//5 0 -17//3; 0 0 17//3 0; 17//3 -17//3 0 0])))
    @test is_isomorphic_form(f1LD2, GWClass(matrix(QQ, [-3//5 0; 0 -17//3])))
end

@testset "M2 Test 6 — quartic over GF(17), local-global" begin
    k = GF(17)
    T3, (w,) = polynomial_ring(k, ["w"])
    f2 = [w^4 + w^3 - w^2 - w]
    f2GD = global_A1_degree(f2)
    @test gw_witt_index(f2GD) == 2
    @test form_rank(f2GD) == 4                 # M2 capture
    f2LD1 = local_A1_degree(f2, ideal(T3, [w + 1]))
    @test gw_witt_index(f2LD1) == 1
    f2LD2 = local_A1_degree(f2, ideal(T3, [w - 1]))
    f2LD3 = local_A1_degree(f2, ideal(T3, [w]))
    f2LDsum = gw_direct_sum(gw_direct_sum(f2LD1, f2LD2), f2LD3)
    @test is_isomorphic_form(f2LDsum, f2GD)
    # M2 captures (class-level):
    @test is_isomorphic_form(f2LD1, GWClass(matrix(k, [1 -1; -1 -3])))
    @test is_isomorphic_form(f2LD2, GWClass(matrix(k, [4;;])))
    @test is_isomorphic_form(f2LD3, GWClass(matrix(k, [-1;;])))
end

@testset "M2 Test 9 — getLocalAlgebraBasis" begin
    S, (x, y) = polynomial_ring(QQ, ["x", "y"])
    f = [x^2 + 1 - y, y]
    p = ideal(S, [x^2 + 1, y])
    B = local_algebra_basis(f, p)
    # M2: {1, x}. Oscar's monomial_basis returns the same monomials in a different
    # order ([x, 1]) — compared as a set.
    @test length(B) == 2
    @test Set(B) == Set([one(S), x])
end

@testset "Doc-example cross-check — cubic meets the x-axis (M2 capture)" begin
    Sd, (u, v) = polynomial_ring(QQ, ["u", "v"])
    f3 = [u^3 - u^2 - v, v]
    g = global_A1_degree(f3)
    l1 = local_A1_degree(f3, ideal(Sd, [u - 1, v]))
    l2 = local_A1_degree(f3, ideal(Sd, [u, v]))
    @test is_isomorphic_form(g, gw_direct_sum(l1, l2))
    @test is_isomorphic_form(g, GWClass(matrix(QQ, [0 -1 1; -1 1 0; 1 0 0])))
    @test is_isomorphic_form(l1, GWClass(matrix(QQ, [1;;])))
    @test is_isomorphic_form(l2, GWClass(matrix(QQ, [0 -1; -1 1])))
    @test form_signature(g) == 1               # doc: signed count of crossings
end

@testset "Extra local-global cross-checks (independent of M2)" begin
    # x^3 - x over QQ: three rational simple zeros
    Sx, (t,) = polynomial_ring(QQ, ["t"])
    f4 = [t^3 - t]
    g4 = global_A1_degree(f4)
    locals4 = [local_A1_degree(f4, ideal(Sx, [t])),
               local_A1_degree(f4, ideal(Sx, [t - 1])),
               local_A1_degree(f4, ideal(Sx, [t + 1]))]
    @test is_isomorphic_form(g4, reduce(gw_direct_sum, locals4))
    @test form_rank(g4) == 3

    # w^2 - 1 over GF(7): two simple zeros
    k7 = GF(7)
    S7, (s,) = polynomial_ring(k7, ["s"])
    f5 = [s^2 - 1]
    g5 = global_A1_degree(f5)
    locals5 = [local_A1_degree(f5, ideal(S7, [s - 1])),
               local_A1_degree(f5, ideal(S7, [s + 1]))]
    @test is_isomorphic_form(g5, reduce(gw_direct_sum, locals5))

    # a zero of multiplicity 2 at a non-rational point: x^2 - 2 over QQ
    f6 = [t^2 - 2]
    g6 = global_A1_degree(f6)
    l6 = local_A1_degree(f6, ideal(Sx, [t^2 - 2]))
    @test is_isomorphic_form(g6, l6)           # single point carries the whole degree
end

@testset "Error paths (M2 LocalGlobalDegrees.m2 source)" begin
    S2, (a, b) = polynomial_ring(QQ, ["a", "b"])
    # non-isolated zeros
    @test_throws ErrorException global_A1_degree([a * b, a * b])
    # number of polynomials != number of variables (zero-dimensional input)
    S1, (t,) = polynomial_ring(QQ, ["t"])
    @test_throws ErrorException global_A1_degree([t^2, t^3])
    # local degree at an ideal the polynomials do not vanish on
    @test_throws ErrorException local_A1_degree([t^3 - t], ideal(S1, [t - 5]))
end
