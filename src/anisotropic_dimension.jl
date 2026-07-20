# Port of M2 Code/AnisotropicDimension.m2: isHyperbolicQQp,
# getAnisotropicDimensionQQp, getAnisotropicDimensionQQ, getAnisotropicDimension,
# getWittIndex. Koprowski-Czogala Algorithms 6, 8, 9.
#
# Note the SIGNED discriminants: Koprowski-Czogala use (-1)^(r(r±1)/2) * disc —
# port the exponents exactly.

const _UNSUPPORTED_FIELD_ERR = "Base field not supported; only implemented over " *
    "QQ, RR, CC, and finite fields of characteristic not 2"

# M2: isHyperbolicQQp — is beta totally hyperbolic over QQ_p? (Algorithm 6)
function _is_hyperbolic_qqp(beta::GWClass, p)
    gw_base_field(beta) == QQ || error("GrothendieckWittClass is not over QQ")
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    r = form_rank(beta)
    isodd(r) && return false
    d = (-1)^div(r * (r - 1), 2) * integral_discriminant(beta)
    _is_padic_square(d, p) || return false
    m = div(r, 2)
    hw_hyperbolic = hilbert_symbol_padic(-1, -1, p)^div(m * (m - 1), 2)
    hw_hyperbolic == hasse_witt_invariant(beta, p)
end

# M2: getAnisotropicDimensionQQp — always 0..4 (Algorithm 8)
"""
    anisotropic_dimension_qqp(beta::GWClass, p)

The anisotropic dimension of a rational form over the ``p``-adic completion
``\\mathbb{Q}_p``: the rank of the anisotropic part of `beta` base-changed to ``\\mathbb{Q}_p``.
Every form of rank ≥ 5 over ``\\mathbb{Q}_p`` is isotropic, so the result is always
0, 1, 2, 3, or 4. This implements [KC18, Algorithm 8].

# Examples
```julia-repl
julia> anisotropic_dimension_qqp(diagonal_form(QQ, (1, -1, 2)), 2)
1
```

# References
- [KC18] P. Koprowski and A. Czogała, *Computing with quadratic forms over number fields*, Journal of Symbolic Computation, 2018.

See also [`anisotropic_dimension`](@ref).
"""
function anisotropic_dimension_qqp(beta::GWClass, p)
    gw_base_field(beta) == QQ || error("GrothendieckWittClass is not over QQ")
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    r = form_rank(beta)
    if iseven(r)
        _is_hyperbolic_qqp(beta, p) && return 0
        d = (-1)^div(r * (r - 1), 2) * integral_discriminant(beta)
        return _is_padic_square(d, p) ? 4 : 2
    else
        c = (-1)^div(r * (r + 1), 2) * integral_discriminant(beta)
        gamma = gw_direct_sum(beta, diagonal_form(QQ, (c,)))
        return _is_hyperbolic_qqp(gamma, p) ? 1 : 3
    end
end

# M2: getAnisotropicDimensionQQ — max over RR (|signature|), QQ_2, and the
# relevant primes (Algorithm 9)
function _anisotropic_dimension_qq(beta::GWClass)
    gw_base_field(beta) == QQ || error("GrothendieckWittClass is not over QQ")
    dims = Int[abs(form_signature(beta)), anisotropic_dimension_qqp(beta, 2)]
    for p in relevant_primes(beta)
        push!(dims, anisotropic_dimension_qqp(beta, p))
    end
    maximum(dims)
end

# M2: getAnisotropicDimension, per field
"""
    anisotropic_dimension(beta::GWClass)
    anisotropic_dimension(A)

The anisotropic dimension of a symmetric bilinear form over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a
finite field of odd characteristic. By the Witt decomposition theorem any
nondegenerate form decomposes uniquely as ``β ≅ n·\\mathbb{H} ⊕ β_a`` with ``β_a``
anisotropic; the anisotropic dimension is the rank of ``β_a``.

Over ``\\mathbb{Q}`` it is the maximum of the anisotropic dimensions over all completions:
``\\lvert \\text{signature} \\rvert`` at the real place and
[`anisotropic_dimension_qqp`](@ref) at 2 and the [`relevant_primes`](@ref)
([KC18, Algorithm 9]).

# Examples
```julia-repl
julia> anisotropic_dimension(diagonal_form(QQ, (1, -1, 2)))
1
```

# References
- [KC18] P. Koprowski and A. Czogała, *Computing with quadratic forms over number fields*, Journal of Symbolic Computation, 2018.

See also [`gw_witt_index`](@ref), [`anisotropic_part`](@ref),
[`is_anisotropic_form`](@ref).
"""
function anisotropic_dimension(A::MatElem)
    R = base_ring(A)
    (R isa QQField || (R isa FinField && characteristic(R) != 2)) ||
        error(_UNSUPPORTED_FIELD_ERR)
    _is_symmetric_square(A) || error("Matrix is not symmetric")
    if R isa QQField
        return _anisotropic_dimension_qq(GWClass(_nondegenerate_part_diagonal(A)))
    else
        diagA = diagonalize_via_congruence(A)
        r = form_rank(diagA)   # number of nonzero diagonal entries
        isodd(r) && return 1
        if is_square(det(_nondegenerate_part_diagonal(diagA))) ==
           is_square(R((-1)^div(r, 2)))
            return 0
        else
            return 2
        end
    end
end

function anisotropic_dimension(A::Matrix{Float64})
    _is_symmetric_square(A) || error("Matrix is not symmetric")
    pos, neg = _pos_neg_diag_counts(diagonalize_via_congruence(A))
    abs(pos - neg)
end

function anisotropic_dimension(A::Matrix{ComplexF64})
    _is_symmetric_square(A) || error("Matrix is not symmetric")
    form_rank(A) % 2
end

anisotropic_dimension(alpha::GWClass) = anisotropic_dimension(gw_matrix(alpha))

# M2: getWittIndex — (rank - anisotropic dimension)/2
"""
    gw_witt_index(beta::GWClass)

The Witt index of a form over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite field of odd
characteristic: the number ``n`` of hyperbolic summands in the Witt
decomposition ``β ≅ n·\\mathbb{H} ⊕ β_a`` ([L05, I.4.3]), computed as
(rank − anisotropic dimension)/2.

The name carries the `gw_` prefix because Oscar exports `witt_index`.

# Examples
```julia-repl
julia> gw_witt_index(diagonal_form(QQ, (1, -1, 2)))
1
```

# References
- [L05] T. Y. Lam, *Introduction to quadratic forms over fields*, American Mathematical Society, 2005.

See also [`anisotropic_dimension`](@ref), [`sum_decomposition`](@ref).
"""
gw_witt_index(alpha::GWClass) =
    div(form_rank(alpha) - anisotropic_dimension(alpha), 2)
