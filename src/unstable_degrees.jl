# Port of M2 Code/UnstableLocalGlobalDegrees.m2 — the exact (QQ, F_q) paths of
# getGlobalUnstableA1Degree / getLocalUnstableA1Degree, plus the internal helpers
# getMultiplicity and makeAntidiagonalUnstableForm (from BuildingForms.m2).
#
# The ℂ paths (roots + linearTolerance cancellation: getLocalUnstableA1DegreeCC,
# the two-input CC variant) are the numerical layer (numerical.jl).
#
# Representation note: rational functions are elements of Oscar's fraction field
# of a one-variable polynomial ring (mirroring M2's frac k[x]); plain polynomials
# are accepted and treated as f/1, as in M2's sub(q, frac R). Oscar fractions
# auto-reduce and auto-normalize the denominator to monic; M2's own normalization
# (f and g both divided by leadCoefficient f) makes the final (f, g) pair
# identical to M2's, since f/g is representative-independent.

const _UNSTABLE_FIELD_ERR =
    "only implemented over QQ and finite fields of characteristic not 2 " *
    "(over CC, pass HomotopyContinuation expressions for the numerical computation)"

# M2: getMultiplicity (internal) — order of vanishing of f at r
function _multiplicity(f::MPolyRingElem, r)
    S = parent(f)
    ngens(S) == 1 || error("need polynomial with one variable")
    u = gens(S)[1]
    mult = 0
    while iszero(evaluate(f, [r]))
        mult += 1
        f = divexact(f, u - r)
    end
    mult
end

# M2: makeAntidiagonalUnstableForm (internal) — n x n antidiagonal with entry b
function _antidiagonal_unstable_form(kk, b, n::Integer)
    A = zero_matrix(kk, n, n)
    for i in 1:n
        A[i, n + 1 - i] = b
    end
    GWuClass(A)
end

# M2: getGlobalUnstableA1Degree (RingElement) — the one-variable Bézoutian of a
# pointed rational function: D = (f(X)g(Y) - f(Y)g(X)) / (X - Y), coefficients
# read against the monomial bases 1, X, ..., X^m and 1, Y, ..., Y^m.
function _global_unstable_A1_degree(f::MPolyRingElem, g::MPolyRingElem)
    S = parent(f)
    kk = coefficient_ring(S)
    (kk isa QQField || (kk isa FinField && characteristic(kk) != 2)) ||
        error(_UNSTABLE_FIELD_ERR)

    # M2 normalizes g by leadCoefficient(f) first, then f
    lc = leading_coefficient(f)
    g = g * inv(lc)
    f = f * inv(lc)

    dim(ideal(S, [f])) > 0 &&
        error("rational function does not have isolated zeros")
    ngens(S) == 1 ||
        error("the number of variables does not match the number of polynomials")
    total_degree(f) > total_degree(g) ||
        error("the rational function is not pointed")

    R2, (X, Y) = polynomial_ring(kk, ["X", "Y"])
    hX = hom(S, R2, [X])
    hY = hom(S, R2, [Y])
    D = divexact(hX(f) * hY(g) - hY(f) * hX(g), X - Y)

    m = degree(D, 1)
    n = degree(D, 2)
    B = zero_matrix(kk, m + 1, m + 1)
    for i in 0:m, j in 0:n
        B[i + 1, j + 1] = coeff(D, [i, j])
    end
    GWuClass(B)
end

"""
    global_unstable_A1_degree(q)
    global_unstable_A1_degree(f, g)

The global unstable ``\\mathbb{A}^1``-Brouwer degree of a pointed rational function
``f/g : \\mathbb{P}^1_k → \\mathbb{P}^1_k`` — pointed meaning ``(f/g)(∞) = ∞``, i.e.
``\\deg f > \\deg g`` — as a [`GWuClass`](@ref) in the unstable
Grothendieck–Witt group
``\\text{GW}^u(k) = \\text{GW}(k) ×_{k^×/(k^×)^2} k^×``.

Morel's ``\\mathbb{A}^1``-Brouwer degree generalizes the classical Brouwer degree by
assigning to an endomorphism of the motivic sphere a class in the
Grothendieck–Witt ring. That degree map is an isomorphism in dimensions two
and above, but in dimension one it is only surjective [M12]; there, a
computation of Morel [M12] and Cazanave [C12] refines it to an isomorphism
``[\\mathbb{P}^1_k, \\mathbb{P}^1_k] ≅ \\text{GW}^u(k)`` onto the unstable group, which records not
only the stable class but also a ``k^×``-scalar. Building on Cazanave's work,
Kass–Wickelgren [KW20] and Igieobo et al. [I+24] give an explicit bilinear
form representing the degree of ``f/g`` in both the local and global settings,
a variant of the Bézoutian form (Cazanave [C12, Thm. 3.6]).

Unlike the stable degree, the global unstable degree is **not** the sum of the
[`local_unstable_A1_degree`](@ref)s at the zeros of ``f/g``: it is their
[`divisorial_sum`](@ref) [I+24], which weights each zero's contribution by the
configuration of the whole divisor of zeros.

# Input and base field

`q` is an element of the fraction field of a one-variable polynomial ring over
``\\mathbb{Q}`` or a finite field of odd characteristic (a plain polynomial is treated
as ``f/1``); the two-argument form supplies numerator and denominator
separately. If ``f`` and ``g`` share a common factor it is cancelled and the
reduced function checked for pointedness before the degree is computed. Over
``\\mathbb{R}``, compute over ``\\mathbb{Q}`` and base-change.

Over ``\\mathbb{C}`` — the one case that needs numerical computation, since the class
carries a ``k^×``-scalar and so the actual complex roots must be found — pass
two `HomotopyContinuation` expressions `f`, `g`. Roots of `f` and `g` closer
than the `linear_tolerance` keyword (default `1e-6`) are treated as a common
factor and cancelled.

# Examples

A degree-5 pointed rational function; its rank equals the number of zeros of
``f/g`` counted with multiplicity over ``\\mathbb{C}``:

```julia-repl
julia> S, (x,) = polynomial_ring(QQ, ["x"]);

julia> q = (x^5 - 6*x^4 + 11*x^3 - 2*x^2 - 12*x + 8) // (x^4 - 5*x^2 + 7*x + 1);

julia> global_unstable_A1_degree(q)
([-68 38 11 -14 1; 38 -63 63 -29 7; 11 63 -84 39 -5; -14 -29 39 -16 0; 1 7 -5 0 1], -53240)
```

The divisorial sum of the local degrees at the zeros ``-1, 1, 2`` recovers the
global degree:

```julia-repl
julia> degs = [local_unstable_A1_degree(q, r) for r in [-1, 1, 2]];

julia> is_isomorphic_form(divisorial_sum(degs, [-1, 1, 2]),
                          global_unstable_A1_degree(q))
true
```

The same computation over ``\\mathbb{C}`` with HomotopyContinuation input:

```julia-repl
julia> import HomotopyContinuation; HomotopyContinuation.@var x;

julia> global_unstable_A1_degree((x-1)*(x-2)*(x-3), (x-1)*(x-4))
(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im], -1.9999999999999987 + 3.337899271044116e-15im)
```

# References

- [M12] F. Morel, *``\\mathbb{A}^1``-algebraic topology over a field*, Springer Lecture Notes in Mathematics, 2012.
- [C12] C. Cazanave, *Algebraic homotopy classes of rational functions*, Annales Scientifiques de l'École Normale Supérieure, 2012.
- [KW20] J. L. Kass and K. Wickelgren, *A classical proof that the algebraic homotopy class of a rational function is the residue pairing*, Linear Algebra and its Applications, 2020.
- [I+24] J. Igieobo et al., *Motivic configurations on the line*, Advances in Mathematics 482 (2025), 110637.

See also [`local_unstable_A1_degree`](@ref), [`divisorial_sum`](@ref),
[`global_A1_degree`](@ref).
"""
global_unstable_A1_degree(q::FracElem{<:MPolyRingElem}) =
    _global_unstable_A1_degree(numerator(q), denominator(q))
global_unstable_A1_degree(q::MPolyRingElem) =
    _global_unstable_A1_degree(q, one(parent(q)))
# M2's two-input variant delegates to f/g over exact fields (CC version: numerical.jl)
global_unstable_A1_degree(f::MPolyRingElem, g::MPolyRingElem) =
    global_unstable_A1_degree(f // g)

# M2: getLocalUnstableA1Degree — at a root r of f/g, the form is the m x m
# antidiagonal with entry F(r), where m is the multiplicity of r and
# F = (u - r)^m g / f.
"""
    local_unstable_A1_degree(q, r)
    local_unstable_A1_degree(f, g, r)

The local unstable ``\\mathbb{A}^1``-Brouwer degree of a pointed rational function
``f/g : \\mathbb{P}^1_k → \\mathbb{P}^1_k`` at a zero `r` in the base field, as a
[`GWuClass`](@ref) in ``\\text{GW}^u(k)``. If `r` is a zero of multiplicity ``m``,
the result is the ``m × m`` antidiagonal form with entry the value of
``(u - r)^m · g/f`` at ``r``.

Input shapes match [`global_unstable_A1_degree`](@ref) (see there for
background and references): a fraction or polynomial plus the root, or
numerator and denominator separately; the numerical ``\\mathbb{C}`` path takes two
`HomotopyContinuation` expressions and a number. Non-reduced input is reduced
(and re-checked for pointedness) first.

# Examples
```julia-repl
julia> S, (x,) = polynomial_ring(QQ, ["x"]);

julia> local_unstable_A1_degree((x^2 + x - 2) // (3*x + 5), -2)
([1//3], 1//3)
```

See also [`global_unstable_A1_degree`](@ref), [`divisorial_sum`](@ref),
[`local_A1_degree`](@ref).
"""
function local_unstable_A1_degree(q::FracElem{<:MPolyRingElem}, r)
    S = base_ring(parent(q))
    kk = coefficient_ring(S)
    (kk isa QQField || (kk isa FinField && characteristic(kk) != 2)) ||
        error(_UNSTABLE_FIELD_ERR)
    rr = _coerce_into(kk, r)
    rr === nothing && error("root not from the base field of the polynomial")
    ngens(S) == 1 || error("must input function of one variable")

    f = numerator(q)
    g = denominator(q)
    dim(ideal(S, [f])) > 0 &&
        error("rational function does not have isolated zeros")
    iszero(evaluate(f, [rr])) ||
        error("the field element is not a zero of the function")
    total_degree(f) > total_degree(g) ||
        error("the rational function is not pointed")

    m = _multiplicity(f, rr)
    u = gens(S)[1]
    F = (u - rr)^m * g // f          # auto-reduced: denominator nonzero at rr
    Fr = evaluate(numerator(F), [rr]) * inv(evaluate(denominator(F), [rr]))
    _antidiagonal_unstable_form(kk, Fr, m)
end

local_unstable_A1_degree(q::MPolyRingElem, r) =
    local_unstable_A1_degree(q // one(parent(q)), r)
# M2's three-input (f, g, r) variant delegates to f/g over exact fields
local_unstable_A1_degree(f::MPolyRingElem, g::MPolyRingElem, r) =
    local_unstable_A1_degree(f // g, r)
