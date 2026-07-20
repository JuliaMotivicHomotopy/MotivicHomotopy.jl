# Core types (GWClass / GWuClass).
#
# Ported from the TEST blocks of the A1BrouwerDegrees Macaulay2 package (v2.0).
# Expected values are M2's — taken from the M2 TEST blocks, not derived independently.
#
# Porting notes:
# - M2 `===` on classes is literal Gram-matrix (+ scalar) equality with the cache
#   ignored (CacheTables are ignored by M2 strict comparison); here that is Julia
#   `==` on our types.
# - M2 RR_53 / CC_53 are represented by Matrix{Float64} / Matrix{ComplexF64};
#   gw_base_field / gw_algebra return Float64 / ComplexF64 for these.
# - Assertions inside the same M2 TEST blocks that need getDiagonalClass,
#   isIsomorphicForm, or degree computations are ported in the other test files;
#   each omission is marked below.

using Test
using Oscar
using MotivicHomotopy

@testset "M2 Test 3 — GWClass basics, addGW, multiplyGW over QQ" begin
    M1 = matrix(QQ, [1 0; 0 1])
    M2 = matrix(QQ, [1 2; 2 5])
    G1 = GWClass(M1)
    G2 = GWClass(M2)
    @test gw_base_field(G1) == QQ
    @test gw_matrix(G1) == M1
    G3 = gw_direct_sum(G1, G2)        # M2: addGW
    G4 = gw_tensor_product(G1, G2)    # M2: multiplyGW
    @test gw_matrix(G3) == matrix(QQ, [1 0 0 0; 0 1 0 0; 0 0 1 2; 0 0 2 5])
    @test gw_matrix(G4) == matrix(QQ, [1 2 0 0; 2 5 0 0; 0 0 1 2; 0 0 2 5])
    # operator aliases (thin wrappers; named functions are canonical)
    @test G1 + G2 == G3
    @test G1 * G2 == G4
end

@testset "M2 Test 18 — makeGWClass over QQ and GF(7)" begin
    M1 = matrix(QQ, [1 0 0; 0 1 0; 0 0 1])
    M2 = matrix(QQ, [1 24//10 0; 24//10 -5 0; 0 0 69])
    M3 = matrix(GF(7), [1 0 0; 0 2 0; 0 0 -3])
    @test GWClass(M1) isa GWClass
    @test GWClass(M2) isa GWClass
    @test GWClass(M3) isa GWClass
end

@testset "M2 Test 19 — getBaseField over QQ, RR, CC, GF(7)" begin
    G1 = GWClass(matrix(QQ, [1 0 0; 0 2 3; 0 3 1]))
    G2 = GWClass([1.0 2.4 -2.41; 2.4 -5 0; -2.41 0 69])
    G3 = GWClass(ComplexF64[1im 2.4 -2.41; 2.4 -5 0; -2.41 0 69+1im])
    G4 = GWClass(matrix(GF(7), [1 0 0; 0 2 0; 0 0 -3]))
    @test gw_base_field(G1) == QQ
    @test gw_base_field(G2) == Float64       # M2: instance(getBaseField M2, RealField)
    @test gw_base_field(G3) == ComplexF64    # M2: instance(getBaseField M3, ComplexField)
    @test order(gw_base_field(G4)) == 7      # M2: (getBaseField M4).order == 7
end

@testset "M2 Test 20 — addGW over QQ, RR, CC (literal equality)" begin
    A1 = GWClass(matrix(QQ, [1 0 -3; 0 23 0; -3 0 -2//5]))
    A2 = GWClass(matrix(QQ, [0 1//2 0; 1//2 5//9 0; 0 0 1]))
    A3 = GWClass(matrix(QQ, [1 0 -3 0 0 0; 0 23 0 0 0 0; -3 0 -2//5 0 0 0;
                             0 0 0 0 1//2 0; 0 0 0 1//2 5//9 0; 0 0 0 0 0 1]))
    @test gw_direct_sum(A1, A2) == A3        # M2: addGW(A1, A2) === A3

    B1 = GWClass([sqrt(2) 0 -3; 0 sqrt(5) 0; -3 0 -1/5])
    B2 = GWClass(fill(1/3, 1, 1))
    B3 = GWClass([sqrt(2) 0 -3 0; 0 sqrt(5) 0 0; -3 0 -1/5 0; 0 0 0 1/3])
    @test gw_direct_sum(B1, B2) == B3

    C1 = GWClass(ComplexF64[2im 0 0; 0 -2 0; 0 0 -3])
    C2 = GWClass(ComplexF64[1 0 -3+1im 0; 0 -2 0 0; -3+1im 0 -3 0; 0 0 0 5])
    C3 = GWClass(ComplexF64[2im 0 0 0 0 0 0; 0 -2 0 0 0 0 0; 0 0 -3 0 0 0 0;
                            0 0 0 1 0 -3+1im 0; 0 0 0 0 -2 0 0;
                            0 0 0 -3+1im 0 -3 0; 0 0 0 0 0 0 5])
    @test gw_direct_sum(C1, C2) == C3
end

@testset "M2 Tests 26–28 (stable-part assertions only) — getGWClass" begin
    # The getDiagonalClass / getScalar-of-diagonal assertions of Tests 26-28 are in
    # test_matrix_layer.jl; only the getGWClass(G) === makeGWClass(M) assertions are
    # ported here.
    MR = [0.0 1; 1 0]                              # Test 26 (RR)
    @test stable_part(GWuClass(MR)) == GWClass(MR)
    MC = ComplexF64[1 2 3; 2 4 5; 3 5 7]           # Test 27 (CC)
    @test stable_part(GWuClass(MC)) == GWClass(MC)
    MQ = matrix(QQ, [1 2 3; 2 4 5; 3 5 7])         # Test 28 (QQ)
    @test stable_part(GWuClass(MQ)) == GWClass(MQ)
end

@testset "M2 Test 29 — GWuClass basics and addGWu over QQ" begin
    M1 = matrix(QQ, [1 0; 0 1])
    M2 = matrix(QQ, [1 2; 2 5])
    G1 = GWuClass(M1)
    G2 = GWuClass(M2)
    @test gw_base_field(G1) == QQ
    @test gw_matrix(G1) == M1
    G3 = gw_direct_sum(G1, G2)        # M2: addGWu
    @test gw_matrix(G3) == matrix(QQ, [1 0 0 0; 0 1 0 0; 0 0 1 2; 0 0 2 5])
    @test gw_scalar(G3) == det(M1) * det(M2)
end

@testset "M2 Test 30 — GWuClass constructors and rejections" begin
    M1 = matrix(QQ, [1 0 0; 0 1 0; 0 0 1])
    M2 = matrix(QQ, [1 24//10 0; 24//10 -5 0; 0 0 69])
    M3 = matrix(GF(7), [1 0 0; 0 2 0; 0 0 -3])
    @test GWuClass(M1) isa GWuClass
    @test GWuClass(M2) isa GWuClass
    @test GWuClass(M3) isa GWuClass
    # M2: assert(try(makeGWuClass(M1,-1)) then false else true)  — wrong square class
    @test_throws ErrorException GWuClass(M1, -1)
    # M2: makeGWuClass(M2, sub(3, GF 17)) — scalar from the wrong field entirely
    @test_throws ErrorException GWuClass(M2, GF(17)(3))
    # M2: makeGWuClass(M3, sub(-6, GF 5)) — finite-field scalar of the wrong order
    @test_throws ErrorException GWuClass(M3, GF(5)(-6))
end

@testset "M2 Test 31 — GWu getBaseField and getAlgebra" begin
    G1 = GWuClass(matrix(QQ, [1 0 0; 0 2 3; 0 3 1]))
    G2 = GWuClass([1.0 2.4 -2.41; 2.4 -5 0; -2.41 0 69])
    G3 = GWuClass(ComplexF64[1im 2.4 -2.41; 2.4 -5 0; -2.41 0 69+1im])
    G4 = GWuClass(matrix(GF(7), [1 0 0; 0 2 0; 0 0 -3]))
    @test gw_base_field(G1) == QQ
    @test gw_base_field(G2) == Float64       # M2: getBaseField(M2) === RR_53
    @test gw_base_field(G3) == ComplexF64    # M2: getBaseField(M3) === CC_53
    @test order(gw_base_field(G4)) == 7
    @test gw_algebra(G1) == QQ
    @test gw_algebra(G2) == Float64
    @test gw_algebra(G3) == ComplexF64
    @test order(gw_algebra(G4)) == 7
end

@testset "M2 Test 32 — addGWu over QQ, RR, CC" begin
    A1 = GWuClass(matrix(QQ, [1 0 -3; 0 23 0; -3 0 -2//5]))
    A2 = GWuClass(matrix(QQ, [0 1//2 0; 1//2 5//9 0; 0 0 1]))
    A3 = GWuClass(matrix(QQ, [1 0 -3 0 0 0; 0 23 0 0 0 0; -3 0 -2//5 0 0 0;
                              0 0 0 0 1//2 0; 0 0 0 1//2 5//9 0; 0 0 0 0 0 1]))
    @test gw_direct_sum(A1, A2) == A3        # M2: addGWu(A1, A2) === A3

    B1 = GWuClass([sqrt(2) 0 -3; 0 sqrt(5) 0; -3 0 -1/5])
    B2 = GWuClass(fill(1/3, 1, 1))
    B3 = GWuClass([sqrt(2) 0 -3 0; 0 sqrt(5) 0 0; -3 0 -1/5 0; 0 0 0 1/3])
    @test gw_matrix(gw_direct_sum(B1, B2)) == gw_matrix(B3)
    # M2 compares these scalars only up to 1e-15 (floating-point det of the block sum):
    @test abs(gw_scalar(gw_direct_sum(B1, B2)) - gw_scalar(B3)) < 1e-15

    C1 = GWuClass(ComplexF64[2im 0 0; 0 -2 0; 0 0 -3])
    C2 = GWuClass(ComplexF64[1 0 -3+1im 0; 0 -2 0 0; -3+1im 0 -3 0; 0 0 0 5])
    C3 = GWuClass(ComplexF64[2im 0 0 0 0 0 0; 0 -2 0 0 0 0 0; 0 0 -3 0 0 0 0;
                             0 0 0 1 0 -3+1im 0; 0 0 0 0 -2 0 0;
                             0 0 0 -3+1im 0 -3 0; 0 0 0 0 0 0 5])
    @test gw_direct_sum(C1, C2) == C3        # M2: addGWu(C1, C2) === C3 (exact)
end

@testset "M2 Test 33 — constructors over finite étale algebras" begin
    Rp, (x,) = polynomial_ring(QQ, ["x"])
    A, _ = quo(Rp, ideal(Rp, [x^2 + 1]))     # M2: R = QQ[x]/(x^2 + 1)
    Sp, (y,) = polynomial_ring(QQ, ["y"])
    S, _ = quo(Sp, ideal(Sp, [y^2 - 1]))     # M2: S = QQ[y]/(y^2 - 1)
    M1 = matrix(A, 2, 2, [A(1), A(2), A(2), A(x)])
    M2 = matrix(S, 2, 2, [S(1), S(2), S(2), S(y)])
    @test GWuClass(M1) isa GWuClass          # M2: try(makeGWuClass M1) then true
    @test GWClass(M1) isa GWClass
    @test GWuClass(M2) isa GWuClass
    @test GWClass(M2) isa GWClass
    @test gw_algebra(GWuClass(M1)) === A     # M2: getAlgebra(makeGWuClass M1) === R
    @test gw_algebra(GWClass(M2)) === S      # M2: getAlgebra(makeGWClass M2) === S
end

@testset "addGWuDivisorial — M2 doc example + error paths" begin
    # From UnstableGrothendieckWittClassesDoc.m2 (addGWuDivisorial has no TEST block;
    # M2 Tests 36/37 exercise it only through degree computations, test_v2_features.jl).
    # Expected value CONFIRMED by an M2 1.25.11 run: matrix {{1/3, 0}, {0, 8/3}},
    # scalar 8.
    alpha = GWuClass(matrix(QQ, [1//3;;]))
    beta  = GWuClass(matrix(QQ, [8//3;;]))
    ds = divisorial_sum([alpha, beta], [-2, 1])
    @test gw_matrix(ds) == matrix(QQ, [1//3 0; 0 8//3])
    @test gw_scalar(ds) == QQ(8)
    # Error paths ported from the addGWuDivisorial source:
    @test_throws ErrorException divisorial_sum([alpha, beta], [-2])      # length mismatch
    @test_throws ErrorException divisorial_sum(typeof(alpha)[], Int[])   # empty sum
end

@testset "Well-definedness rejections (M2 isWellDefinedGW error paths)" begin
    # No M2 TEST block covers these; they exercise the documented constructor checks
    # in GrothendieckWittClasses.m2 (symmetry, nondegeneracy, characteristic != 2).
    @test_throws ErrorException GWClass(matrix(QQ, [1 2; 3 4]))      # not symmetric
    @test_throws ErrorException GWClass(matrix(QQ, [1 1; 1 1]))      # degenerate
    @test_throws ErrorException GWClass(matrix(GF(2), [0 1; 1 0]))   # characteristic 2
end
