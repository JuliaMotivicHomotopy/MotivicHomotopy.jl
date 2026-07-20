# Port of M2 Code/LocalGlobalDegrees.m2 (getGlobalA1Degree, getLocalA1Degree)
# plus getLocalAlgebraBasis from Code/ArithmeticMethods.m2.
#
# The Bézoutian construction is ported move-for-move (Brazelton-McKean-Pauli):
#  1. double the variables: R = kk[X_1..X_n, Y_1..Y_n] (X's first, as in M2);
#  2. divided differences: D[i,j] = (f_i(Y_1..Y_{j-1}, X_j..X_n)
#                                    - f_i(Y_1..Y_j, X_{j+1}..X_n)) / (X_j - Y_j);
#  3. reduce det(D) modulo I(X) + I(Y) (global) or J(X) + J(Y) (local);
#  4. read coefficients against products of the standard monomial bases.
# Soundness of step 4 relies on: (a) Oscar's default ordering is degrevlex, same
# as M2's GRevLex; (b) I(X) and I(Y) live in disjoint variables, so their Gröbner
# bases union and the standard monomials of the sum are exactly the products.
#
# Substitutions: M2 divides in frac(R) and lifts back — here
# `divexact` does the exact polynomial division directly. M2's `%` is
# `normal_form`. Basis ORDER differs from M2's `basis`, so Gram matrices are
# permutation-congruent to M2's (same class); degree results are compared as
# classes, as M2's own TESTs do.
#
# Field scope: exact fields only — QQ and F_q (odd q). ℝ: M2 errors and directs
# the user to compute over QQ and base-change. ℂ: M2's rank-only shortcut needs
# the algebra dimension over inexact CC — see the numerical layer (numerical.jl).

# M2: getLocalAlgebraBasis — monomial basis of the local algebra
# Q_p(f) = R[x]_p / (f), realized as R/(I : sat(I, p)).
"""
    local_algebra_basis(L, p)

A monomial basis of the local algebra ``Q_p(f) = k[x_1,…,x_n]_{\\mathfrak{m}_p}/(f)`` of
an endomorphism of affine space at an isolated zero: `L` is the list of
polynomials ``f = (f_1, …, f_n)`` and `p` the prime ideal of the zero. The
local algebra is realized as ``k[x]/(I : (I : p^∞))`` ([S02, Proposition 2.5])
and its standard monomials are returned.

# Examples
```julia-repl
julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> local_algebra_basis([x^2 + 1 - y, y], ideal(S, [x^2 + 1, y]))
2-element Vector{QQMPolyRingElem}:
 x
 1
```

# References
- [S02] B. Sturmfels, *Solving Systems of Polynomial Equations*, American Mathematical Society, 2002.

See also [`local_A1_degree`](@ref).
"""
function local_algebra_basis(L::Vector{<:MPolyRingElem}, p::MPolyIdeal)
    R = parent(L[1])
    I = ideal(R, L)
    coefficient_ring(R) isa Field ||
        error("this computation is only supported over polynomial rings over fields")
    is_subset(I, p) ||
        error("the polynomials of the list do not vanish at the prescribed ideal")
    try
        is_prime(p)
    catch
        println("warning: Unable to verify whether the prime ideal is an " *
                "isolated zero of the list of polynomials")
    end
    dim(I) > 0 && error("morphism does not have isolated zeros")
    is_subset(I, p) || error("prime is not a zero of function")
    J = quotient(I, saturation(I, p))
    monomial_basis(quo(R, J)[1])
end

# Variable doubling + divided differences + Bézoutian determinant.
function _bezoutian_data(Endo::Vector{<:MPolyRingElem}, S, kk)
    n = length(Endo)
    R, XY = polynomial_ring(kk, vcat(["X" * string(i) for i in 1:n],
                                     ["Y" * string(i) for i in 1:n]))
    Xv = XY[1:n]
    Yv = XY[(n + 1):(2 * n)]
    D = Matrix{elem_type(R)}(undef, n, n)
    for j in 1:n
        phi1 = hom(S, R, [k < j ? Yv[k] : Xv[k] for k in 1:n])
        phi2 = hom(S, R, [k <= j ? Yv[k] : Xv[k] for k in 1:n])
        for i in 1:n
            D[i, j] = divexact(phi1(Endo[i]) - phi2(Endo[i]), Xv[j] - Yv[j])
        end
    end
    (R, Xv, Yv, det(matrix(R, D)))
end

function _degree_base_checks(Endo::Vector{<:MPolyRingElem})
    S = parent(Endo[1])
    kk = coefficient_ring(S)
    (kk isa QQField || (kk isa FinField && characteristic(kk) != 2)) ||
        error(_DEGREE_FIELD_ERR)
    (S, kk)
end

const _DEGREE_FIELD_ERR = "Base field not supported; pass exact polynomials " *
    "over QQ or a finite field of odd characteristic, floating-point input " *
    "for the numerical computation over CC, or — over RR — compute over QQ " *
    "and base-change the result"

# Read the Gram matrix of the Bézoutian form: coefficients of the reduced
# determinant against products of the X- and Y-standard-monomial bases.
function _bezoutian_gram(kk, R, Xv, Yv, bez_red, basisX, basisY, RX, RY)
    hXR = hom(RX, R, Xv)
    hYR = hom(RY, R, Yv)
    bx = [hXR(b) for b in basisX]
    by = [hYR(b) for b in basisY]
    m = length(bx)
    B = zero_matrix(kk, m, m)
    for i in 1:m, j in 1:m
        B[i, j] = coeff(bez_red, bx[i] * by[j])
    end
    B
end

# M2: getGlobalA1Degree
"""
    global_A1_degree(F)

The global ``\\mathbb{A}^1``-Brouwer degree of an endomorphism of affine space
``f = (f_1, …, f_n) : \\mathbb{A}^n_k → \\mathbb{A}^n_k`` with isolated zeros, as a
[`GWClass`](@ref) in ``\\text{GW}(k)``. `F` is a vector of ``n`` polynomials
in ``n`` variables over a field ``k`` of characteristic not 2.

The ``\\mathbb{A}^1``-Brouwer degree, first defined by Morel [M12], is an
algebro-geometric enrichment of the classical topological Brouwer degree.
Using the tools of motivic homotopy theory one associates to an endomorphism
of affine space the isomorphism class of a nondegenerate symmetric bilinear
form whose invariants encode geometric data about how the morphism transforms
space: its rank recovers the degree of the associated complex map, and its
signature the degree of the associated real map. Such a form appears in the
work of Eisenbud–Levine [EL77] and Khimshiashvili [K77], whose signature
computes the local degree of a smooth map of real manifolds even where the
Jacobian vanishes; this was shown to agree with Morel's degree by
Kass–Wickelgren [KW19]. A related form attached to a complete intersection,
due to Scheja–Storch [SS76], was aligned with the ``\\mathbb{A}^1``-degree in [BW23].
Following Brazelton–McKean–Pauli [BMP23], the degree is computed here as a
multivariate Bézoutian bilinear form.

Following McKean [M21] one may read the degree ``\\deg^{\\mathbb{A}^1}(f)`` as a
quadratically enriched intersection multiplicity of the hypersurfaces
``V(f_1) ∩ ⋯ ∩ V(f_n)``. It equals the sum of the [`local_A1_degree`](@ref)s
over the points of the zero locus ``V(f)``.

The Gram matrix is expressed in a standard-monomial basis; an equivalent form
in a different basis represents the same class, so compare results with
[`is_isomorphic_form`](@ref) rather than entrywise.

# Base field

- **``\\mathbb{Q}`` and finite fields of odd characteristic** — the form is computed
  exactly.
- **``\\mathbb{C}``** — a symmetric bilinear form over ``\\mathbb{C}`` is determined by its rank,
  so the degree is simply the identity form of rank equal to the
  ``\\mathbb{C}``-dimension of the coordinate algebra
  ``k[x_1,…,x_n]/(f_1,…,f_n)`` — the number of zeros counted with
  multiplicity. This dimension is a discrete invariant equal to the rank of
  the degree computed over ``\\mathbb{Q}``, so it is obtained from the exact
  computation with no numerical root-finding.
- **``\\mathbb{R}``** — compute the degree over ``\\mathbb{Q}`` and base-change the resulting
  Gram matrix to ``\\mathbb{R}`` (its [`form_signature`](@ref) is the real degree);
  real input is not accepted directly.

# Examples

For ``z ↦ z^2`` the degree is a rank-2 form of signature 0: the complex map
``\\mathbb{C}`` → ``\\mathbb{C}`` has degree 2, while the real map ``\\mathbb{R}`` → ``\\mathbb{R}`` has degree 0.

```julia-repl
julia> S, (x,) = polynomial_ring(QQ, ["x"]);

julia> beta = global_A1_degree([x^2 + 1])
[0   1]
[1   0]

julia> form_rank(beta), form_signature(beta)
(2, 0)
```

The cubic ``y = x(x-1)(x+1)`` meeting the ``x``-axis, read as an enriched
count of intersection points: rank 3 (three complex intersections) and
signature 1 (the signed real count).

```julia-repl
julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> f = [x^3 - x^2 - y, y];

julia> global_A1_degree(f)
[0    0    1]
[0    1   -1]
[1   -1    0]

julia> form_signature(global_A1_degree(f))
1
```

The global degree is the sum of the local degrees over the zero locus
``V(f) = \\{(1,0), (0,0)\\}``:

```julia-repl
julia> d1 = local_A1_degree(f, ideal(S, [x - 1, y]));

julia> d2 = local_A1_degree(f, ideal(S, [x, y]));

julia> is_isomorphic_form(global_A1_degree(f), gw_direct_sum(d1, d2))
true
```

# References

- [M12] F. Morel, *``\\mathbb{A}^1``-algebraic topology over a field*, Springer Lecture Notes in Mathematics, 2012.
- [EL77] D. Eisenbud and H. Levine, *An algebraic formula for the degree of a C∞ map germ*, Annals of Mathematics, 1977.
- [K77] G. Khimshiashvili, *The local degree of a smooth mapping*, Sakharth. SSR Mecn. Akad. Moambe, 1977.
- [SS76] G. Scheja and U. Storch, *Über Spurfunktionen bei vollständigen Durchschnitten*, J. Reine Angew. Math., 1975.
- [KW19] J. L. Kass and K. Wickelgren, *The class of Eisenbud–Khimshiashvili–Levine is the local ``\\mathbb{A}^1``-Brouwer degree*, Duke Mathematical Journal, 2019.
- [BW23] T. Bachmann and K. Wickelgren, *Euler classes: six-functors formalism, dualities, integrality and linear subspaces of complete intersections*, J. Inst. Math. Jussieu, 2023.
- [BMP23] T. Brazelton, S. McKean, and S. Pauli, *Bézoutians and the ``\\mathbb{A}^1``-degree*, Algebra & Number Theory, 2023.
- [M21] S. McKean, *An arithmetic enrichment of Bézout's Theorem*, Mathematische Annalen, 2021.

See also [`local_A1_degree`](@ref), [`global_unstable_A1_degree`](@ref),
[`sum_decomposition`](@ref).
"""
function global_A1_degree(Endo::Vector{<:MPolyRingElem})
    n = length(Endo)
    S, kk = _degree_base_checks(Endo)
    dim(ideal(S, Endo)) > 0 && error("morphism does not have isolated zeros")
    ngens(S) == n ||
        error("the number of variables does not match the number of polynomials")

    R, Xv, Yv, bez = _bezoutian_data(Endo, S, kk)

    RX, Xs = polynomial_ring(kk, ["X" * string(i) for i in 1:n])
    RY, Ys = polynomial_ring(kk, ["Y" * string(i) for i in 1:n])
    phiX = hom(S, RX, Xs)
    phiY = hom(S, RY, Ys)
    basisX = monomial_basis(quo(RX, ideal(RX, [phiX(f) for f in Endo]))[1])
    basisY = monomial_basis(quo(RY, ideal(RY, [phiY(f) for f in Endo]))[1])

    hSX = hom(S, R, Xv)
    hSY = hom(S, R, Yv)
    promoted = ideal(R, vcat([hSX(f) for f in Endo], [hSY(f) for f in Endo]))
    bez_red = normal_form(bez, promoted)

    GWClass(_bezoutian_gram(kk, R, Xv, Yv, bez_red, basisX, basisY, RX, RY))
end

"""
    local_A1_degree(F, p)

The local ``\\mathbb{A}^1``-Brouwer degree of an endomorphism of affine space
``f = (f_1, …, f_n) : \\mathbb{A}^n_k → \\mathbb{A}^n_k`` at an isolated zero, as a
[`GWClass`](@ref) in ``\\text{GW}(k)``. `F` is a vector of ``n`` polynomials
in ``n`` variables over a field of characteristic not 2, and `p` is the prime
ideal of a point in the zero locus ``V(f)``.

The local degree is the class of the Bézoutian bilinear form on the local
algebra ``Q_p(f) = k[x_1,…,x_n]_{\\mathfrak{m}_p}/(f_1,…,f_n)`` at the point (see
[`local_algebra_basis`](@ref)). Summed over the points of ``V(f)`` it
recovers the [`global_A1_degree`](@ref); see there for the background and
references.

The base field is handled as for the global degree: exactly over ``\\mathbb{Q}`` and
finite fields of odd characteristic; over ``\\mathbb{C}`` the class is the identity
form of rank equal to the ``\\mathbb{C}``-dimension of ``Q_p(f)`` (the multiplicity of
the zero), obtained from the exact computation; over ``\\mathbb{R}``, compute over
``\\mathbb{Q}`` and base-change.

# Examples
```julia-repl
julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> f = [x^3 - x^2 - y, y];

julia> d1 = local_A1_degree(f, ideal(S, [x - 1, y]))
[1]

julia> d2 = local_A1_degree(f, ideal(S, [x, y]))
[ 1   -1]
[-1    0]

julia> is_isomorphic_form(global_A1_degree(f), gw_direct_sum(d1, d2))
true
```

See also [`global_A1_degree`](@ref), [`local_unstable_A1_degree`](@ref),
[`local_algebra_basis`](@ref).
"""
function local_A1_degree(Endo::Vector{<:MPolyRingElem}, p::MPolyIdeal)
    n = length(Endo)
    S, kk = _degree_base_checks(Endo)
    I = ideal(S, Endo)
    # k[x_1..x_n]_p / I_p is isomorphic to k[x_1..x_n] / J
    J = quotient(I, saturation(I, p))
    dim(I) > 0 && error("morphism does not have isolated zeros")
    ngens(S) == n ||
        error("the number of variables does not match the number of polynomials")

    R, Xv, Yv, bez = _bezoutian_data(Endo, S, kk)

    RX, Xs = polynomial_ring(kk, ["X" * string(i) for i in 1:n])
    RY, Ys = polynomial_ring(kk, ["Y" * string(i) for i in 1:n])
    phiX = hom(S, RX, Xs)
    phiY = hom(S, RY, Ys)
    basisX = local_algebra_basis([phiX(f) for f in Endo],
                                 ideal(RX, [phiX(g) for g in gens(p)]))
    basisY = local_algebra_basis([phiY(f) for f in Endo],
                                 ideal(RY, [phiY(g) for g in gens(p)]))

    hSX = hom(S, R, Xv)
    hSY = hom(S, R, Yv)
    local_ideal = ideal(R, vcat([hSX(g) for g in gens(J)], [hSY(g) for g in gens(J)]))
    bez_red = normal_form(bez, local_ideal)

    GWClass(_bezoutian_gram(kk, R, Xv, Yv, bez_red, basisX, basisY, RX, RY))
end
