# Port of M2 Code/IsomorphismOfForms.m2: isIsomorphicFormQQ (internal) and
# isIsomorphicForm for matrices, GW classes, and unstable GW classes.
#
# Per-field classification (M2's logic):
#   CC: same dimension + same rank
#   RR: same dimension + same positive/negative diagonal counts
#   QQ: same dimension + rank/signature/discriminant/Hasse-Witt of the
#       nondegenerate parts (Hasse-Minkowski)
#   GF(q), q odd: same dimension + same rank + determinants in the same square class
#
# linear_tolerance is used ONLY by the unstable RR/CC scalar comparison — it never
# touches the exact QQ / GF paths. M2 validates it only in the unstable variant;
# ported as-is.

const _DIFFERENT_FIELDS_ERR = "Base fields are not the same"
const _NOT_SYMMETRIC_ERR = "Underlying matrix is not symmetric"

# M2: isIsomorphicFormQQ
function _is_isomorphic_form_qq(alpha::GWClass, beta::GWClass)
    gw_base_field(alpha) == QQ || error("first input must have base field QQ")
    gw_base_field(beta) == QQ || error("second input must have base field QQ")
    form_rank(alpha) == form_rank(beta) || return false
    form_signature(alpha) == form_signature(beta) || return false
    integral_discriminant(alpha) == integral_discriminant(beta) || return false
    for p in unique(vcat(relevant_primes(alpha), relevant_primes(beta)))
        hasse_witt_invariant(alpha, p) == hasse_witt_invariant(beta, p) ||
            return false
    end
    true
end

_supported_exact_field(R) = R isa QQField || (R isa FinField && characteristic(R) != 2)

# transfer x into k1 when both finite fields have the same order (prime fields via ZZ)
function _same_order_transfer(k1::FinField, x::FinFieldElem)
    parent(x) === k1 && return x
    try
        k1(lift(ZZ, x))
    catch
        error("transfer between distinct finite fields of the same order is not " *
              "yet supported")
    end
end

"""
    is_isomorphic_form(alpha, beta; linear_tolerance = 1e-6)

Whether two Grothendieck–Witt classes (or unstable classes, or raw symmetric
matrices) over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite field of odd characteristic represent
the same element of ``\\text{GW}(k)`` (resp. ``\\text{GW}^u(k)``). This is the mathematical
notion of equality; `==` on classes compares Gram matrices literally.

The classification used per field:
- **``\\mathbb{C}``** (and any quadratically closed field): rank alone, since every nonzero
  element is a square.
- **``\\mathbb{R}``**: rank and signature (Sylvester's law of inertia).
- **``\\mathbb{Q}``**: rank, signature, discriminant, and the Hasse–Witt invariants at all
  [`relevant_primes`](@ref) — by the Hasse–Minkowski principle forms over ``\\mathbb{Q}``
  are isomorphic iff they are isomorphic over every completion ([S73, IV
  Thm. 7]; [L05, VI.3.3]). Each Hasse–Witt invariant is a product of values of
  a *symbol* on the diagonal entries ([MH73, III.5.4]).
- **finite fields**: rank and the square class of the discriminant.

For [`GWuClass`](@ref)es the fibered-product structure of ``\\text{GW}^u(k)`` reduces
the test to: stable parts isomorphic and ``k^×``-factors equal. Over ``\\mathbb{Q}`` and
finite fields the scalars must agree exactly; over ``\\mathbb{R}`` and ``\\mathbb{C}`` they are
considered equal when the absolute value of their difference is below
`linear_tolerance` (default `1e-6`).

# Examples
```julia-repl
julia> alpha = GWClass(ComplexF64[2 3 1; 3 -1 0; 1 0 0]);

julia> beta = GWClass(ComplexF64[2 4 -1; 4 5 7; -1 7 9]);

julia> is_isomorphic_form(alpha, beta)
true

julia> is_isomorphic_form(GWClass(QQ[1 4 7; 4 3 2; 7 2 -1]),
                          GWClass(QQ[0 0 1; 0 2 7; 1 7 3]))
false

julia> u1 = GWuClass(QQ[2 3 1; 3 -1 0; 1 0 0], 1);

julia> u2 = GWuClass(QQ[2 3 1; 3 -1 0; 1 0 0], 4);

julia> is_isomorphic_form(u1, u2)    # same stable part, scalars 1 ≠ 4 in ℚ×
false
```

# References
- [S73] J. P. Serre, *A course in arithmetic*, Springer-Verlag, 1973.
- [L05] T. Y. Lam, *Introduction to quadratic forms over fields*, American Mathematical Society, 2005.
- [MH73] J. Milnor and D. Husemoller, *Symmetric bilinear forms*, Springer-Verlag, 1973.

See also [`form_rank`](@ref), [`form_signature`](@ref),
[`integral_discriminant`](@ref), [`hasse_witt_invariant`](@ref).
"""
function is_isomorphic_form(A::MatElem, B::MatElem; linear_tolerance::Real = 1e-6)
    k1, k2 = base_ring(A), base_ring(B)
    _supported_exact_field(k1) || error(_UNSUPPORTED_FIELD_ERR)
    _supported_exact_field(k2) || error(_UNSUPPORTED_FIELD_ERR)
    _is_symmetric_square(A) || error(_NOT_SYMMETRIC_ERR)
    _is_symmetric_square(B) || error(_NOT_SYMMETRIC_ERR)
    if k1 isa QQField && k2 isa QQField
        return nrows(A) == nrows(B) &&
            _is_isomorphic_form_qq(GWClass(_nondegenerate_part_diagonal(A)),
                                   GWClass(_nondegenerate_part_diagonal(B)))
    elseif k1 isa FinField && k2 isa FinField && order(k1) == order(k2)
        nrows(A) == nrows(B) || return false
        form_rank(A) == form_rank(B) || return false
        detB = _same_order_transfer(k1, det(_nondegenerate_part_diagonal(B)))
        return is_square(det(_nondegenerate_part_diagonal(A))) == is_square(detB)
    end
    error(_DIFFERENT_FIELDS_ERR)
end

function is_isomorphic_form(A::Matrix{Float64}, B::Matrix{Float64};
                            linear_tolerance::Real = 1e-6)
    _is_symmetric_square(A) || error(_NOT_SYMMETRIC_ERR)
    _is_symmetric_square(B) || error(_NOT_SYMMETRIC_ERR)
    size(A, 1) == size(B, 1) || return false
    _pos_neg_diag_counts(diagonalize_via_congruence(A)) ==
        _pos_neg_diag_counts(diagonalize_via_congruence(B))
end

function is_isomorphic_form(A::Matrix{ComplexF64}, B::Matrix{ComplexF64};
                            linear_tolerance::Real = 1e-6)
    _is_symmetric_square(A) || error(_NOT_SYMMETRIC_ERR)
    _is_symmetric_square(B) || error(_NOT_SYMMETRIC_ERR)
    size(A, 1) == size(B, 1) && form_rank(A) == form_rank(B)
end

# mixed backends = mixed base fields
is_isomorphic_form(::Union{MatElem, FloatGram}, ::Union{MatElem, FloatGram};
                   linear_tolerance::Real = 1e-6) = error(_DIFFERENT_FIELDS_ERR)

is_isomorphic_form(alpha::GWClass, beta::GWClass; linear_tolerance::Real = 1e-6) =
    is_isomorphic_form(gw_matrix(alpha), gw_matrix(beta))

# M2: isIsomorphicForm for unstable classes — stable parts isomorphic AND the
# kx-factors equal (exactly over QQ/GF; within linear_tolerance over RR/CC).
function is_isomorphic_form(alpha::GWuClass, beta::GWuClass;
                            linear_tolerance::Real = 1e-6)
    linear_tolerance > 0 || error("linearTolerance must be a positive number")
    r1, r2 = gw_scalar(alpha), gw_scalar(beta)
    _scalar_field_supported(r1) || error(_UNSUPPORTED_FIELD_ERR)
    _scalar_field_supported(r2) || error(_UNSUPPORTED_FIELD_ERR)
    if (r1 isa Float64 && r2 isa Float64) || (r1 isa ComplexF64 && r2 isa ComplexF64)
        return abs(r1 - r2) < linear_tolerance &&
            is_isomorphic_form(gw_matrix(alpha), gw_matrix(beta))
    elseif r1 isa QQFieldElem && r2 isa QQFieldElem
        return is_isomorphic_form(gw_matrix(alpha), gw_matrix(beta)) && r1 == r2
    elseif r1 isa FinFieldElem && r2 isa FinFieldElem &&
           order(parent(r1)) == order(parent(r2))
        return is_isomorphic_form(gw_matrix(alpha), gw_matrix(beta)) &&
            r1 == _same_order_transfer(parent(r1), r2)
    end
    error(_DIFFERENT_FIELDS_ERR)
end

_scalar_field_supported(r::Union{Float64, ComplexF64, QQFieldElem}) = true
_scalar_field_supported(r::FinFieldElem) = characteristic(parent(r)) != 2
_scalar_field_supported(::Any) = false
