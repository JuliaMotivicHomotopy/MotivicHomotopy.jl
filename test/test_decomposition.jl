# Decomposition layer: anisotropic part, sum decomposition.
#
# Provenance of expected values:
# - "M2 TEST 47": the only upstream TEST block for this layer (trivial paths only).
# - "M2 capture": captured by running M2 1.25.11 on the package —
#   getSumDecomposition has NO upstream TEST block, and the QQ dim-4/3/2 reducers
#   are only hit via trivial shortcuts upstream, so these runs are the oracle.
# - "consistency": validated via invariants and the Witt-decomposition identity
#   q ≅ a ⊕ wH, not against captured M2 output.
#
# CRT/GF(2)-solver-dependent paths are compared up to isomorphism (different valid
# solver outputs give different — but isomorphic — representatives); solver-free
# paths are compared literally.

using Test
using Oscar
using MotivicHomotopy

# q ≅ (anisotropic part) ⊕ (Witt index copies of H), part anisotropic of full rank
function check_witt_decomposition(q)
    a = anisotropic_part(q)
    @test is_anisotropic_form(a) || form_rank(a) == 0
    @test form_rank(a) == anisotropic_dimension(q)
    @test is_isomorphic_form(q, gw_direct_sum(a, hyperbolic_form(QQ, 2 * gw_witt_index(q))))
    a
end

@testset "M2 Test 47 — getAnisotropicPart (trivial paths)" begin
    H = hyperbolic_form(QQ)
    @test form_rank(anisotropic_part(H)) == 0
    twoPos = diagonal_form(QQ, (1, 1))
    @test is_isomorphic_form(anisotropic_part(twoPos), twoPos)
    twoH = diagonal_form(QQ, (1, -1, 1, -1))
    @test form_rank(anisotropic_part(twoH)) == 0
end

@testset "QQ dim-4 reducer, positive signature (M2 capture)" begin
    q4 = diagonal_form(QQ, (1, 1, 1, 1, 1, -1))
    @test anisotropic_dimension(q4) == 4          # M2: 4
    a = check_witt_decomposition(q4)
    # M2 capture: identity 4x4
    @test is_isomorphic_form(a, GWClass(matrix(QQ, [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])))
end

@testset "QQ dim-4 by p-adic obstruction, full 4→3→2 cascade (M2 capture)" begin
    q4p = diagonal_form(QQ, (1, 1, -3, -3, 1, -1))
    @test anisotropic_dimension(q4p) == 4         # M2: 4 (signature is 0)
    a = check_witt_decomposition(q4p)
    # M2 capture: diag(1, 1, -3, -3)
    @test is_isomorphic_form(a, diagonal_form(QQ, (1, 1, -3, -3)))
end

@testset "QQ dim-3 entry, CRT reducer + dim-2 solver (M2 capture)" begin
    q3 = diagonal_form(QQ, (1, 1, 1, 1, -1))
    @test anisotropic_dimension(q3) == 3          # M2: 3
    a = check_witt_decomposition(q3)
    # M2 capture: identity 3x3
    @test is_isomorphic_form(a, diagonal_form(QQ, (1, 1, 1)))
end

@testset "QQ dim-2 entry with rank > 2, GF(2) solver (M2 capture)" begin
    q2 = diagonal_form(QQ, (3, 5, 1, -1))
    @test anisotropic_dimension(q2) == 2          # M2: 2
    a = check_witt_decomposition(q2)
    # M2 capture: diag(3, 5)
    @test is_isomorphic_form(a, diagonal_form(QQ, (3, 5)))
end

@testset "QQ dim-1 branch (M2 capture, solver-free: literal)" begin
    q1 = diagonal_form(QQ, (2, 1, -1))
    @test anisotropic_dimension(q1) == 1          # M2: 1
    a = check_witt_decomposition(q1)
    @test gw_matrix(a) == matrix(QQ, [2;;])       # M2 capture: {{2}}
end

@testset "Negative-signature reducer paths (consistency)" begin
    # Negative-signature forms whose reduction passes through the dim-4 and dim-3
    # reducers; validated via invariants and the Witt decomposition q ≅ a ⊕ wH.
    q4neg = diagonal_form(QQ, (-1, -1, -1, -1, -1, 1))
    @test anisotropic_dimension(q4neg) == 4
    a = check_witt_decomposition(q4neg)
    @test is_isomorphic_form(a, diagonal_form(QQ, (-1, -1, -1, -1)))

    q3neg = diagonal_form(QQ, (-1, -1, -1, -1, 1))
    @test anisotropic_dimension(q3neg) == 3
    b = check_witt_decomposition(q3neg)
    @test is_isomorphic_form(b, diagonal_form(QQ, (-1, -1, -1)))
end

@testset "Per-field anisotropic parts (M2 capture, literal)" begin
    # RR: diag(1,1,-1) -> [[1]]
    @test gw_matrix(anisotropic_part(GWClass([1.0 0 0; 0 1 0; 0 0 -1]))) == [1.0;;]
    # CC: identity 3 -> [[1]]
    @test gw_matrix(anisotropic_part(GWClass(ComplexF64[1 0 0; 0 1 0; 0 0 1]))) ==
        ComplexF64[1;;]
    # GF(7): <1,1,3> -> [[-3]]; <1,3> -> empty
    @test gw_matrix(anisotropic_part(GWClass(matrix(GF(7), [1 0 0; 0 1 0; 0 0 3])))) ==
        matrix(GF(7), [-3;;])
    @test form_rank(anisotropic_part(GWClass(matrix(GF(7), [1 0; 0 3])))) == 0
end

@testset "getSumDecomposition / getSumDecompositionString (M2 capture)" begin
    # QQ doc example (dim-1 path, deterministic)
    gammaQ = GWClass(matrix(QQ, [1 2 3; 2 4 5; 3 5 6]))
    @test gw_matrix(sum_decomposition(gammaQ)) == matrix(QQ, [1 0 0; 0 1 0; 0 0 -1])
    @test sum_decomposition_string(gammaQ) == "H + <1>"
    @test is_isomorphic_form(sum_decomposition(gammaQ), gammaQ)

    # GF(13) doc example — note M2 prints fp elements with BALANCED representatives
    deltaGF = GWClass(matrix(GF(13), [9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5]))
    @test gw_matrix(sum_decomposition(deltaGF)) ==
        matrix(GF(13), [1 0 0 0; 0 -5 0 0; 0 0 1 0; 0 0 0 -1])
    @test sum_decomposition_string(deltaGF) == "H + <1> + <-5>"

    # CC / RR doc examples
    alphaCC = GWClass(ComplexF64[1 2 3; 2 4 5; 3 5 6])
    @test sum_decomposition_string(alphaCC) == "H + <1>"
    betaRR = GWClass([2.091 2.728 6.747; 2.728 7.329 6.257; 6.747 6.257 0.294])
    @test sum_decomposition_string(betaRR) == "H + <1>"

    # rank 0
    @test sum_decomposition_string(GWClass(zero_matrix(QQ, 0, 0))) == "empty form"
end

@testset "Unstable sum decomposition (M2 capture)" begin
    u = GWuClass(matrix(QQ, [0 1; 1 0]), -4)
    su = sum_decomposition(u)
    @test gw_matrix(su) == matrix(QQ, [1 0; 0 -1])
    @test gw_scalar(su) == QQ(-4)
    @test sum_decomposition_string(u) == "(H, -4)"
end

@testset "Cache overwrite mirrors M2" begin
    # M2 capture: getDiagonalClass before = diag(2,-2); after getSumDecomposition,
    # the SAME call returns the sum decomposition diag(1,-1).
    b = GWClass(matrix(QQ, [0 1; 1 0]))
    @test gw_matrix(diagonal_class(b)) == matrix(QQ, [2 0; 0 -2])
    sum_decomposition(b)
    @test gw_matrix(diagonal_class(b)) == matrix(QQ, [1 0; 0 -1])
    # same for the unstable variant
    bu = GWuClass(matrix(QQ, [0 1; 1 0]))
    sum_decomposition(bu)
    @test gw_matrix(diagonal_class(bu)) == matrix(QQ, [1 0; 0 -1])
end

@testset "Error paths (M2 Decomposition.m2 source)" begin
    @test_throws ErrorException anisotropic_part(matrix(QQ, [1 2; 3 4]))
end
