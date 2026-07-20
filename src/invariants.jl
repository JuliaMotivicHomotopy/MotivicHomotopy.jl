# Port of M2 Code/GWInvariants.m2: getRank, countPos/NegDiagEntries (internal),
# getSignature, getIntegralDiscriminant, getRelevantPrimes, getHasseWittInvariant.

# --- M2: getRank ------------------------------------------------------------------
# On a class: the matrix size. On a raw matrix: the actual rank (drops degeneracy).

"""
    form_rank(beta::GWClass)
    form_rank(M)

The rank of a symmetric bilinear form. On a [`GWClass`](@ref) (which is
nondegenerate by construction) this is the size of the Gram matrix; on a raw
matrix it is the matrix rank, so degenerate directions are not counted.

# Examples
```julia-repl
julia> form_rank(diagonal_form(QQ, (3, 5, 7, 11)))
4
```

See also [`form_signature`](@ref), [`anisotropic_dimension`](@ref).
"""
form_rank(beta::GWClass) = _n_of(beta.matrix)
form_rank(M::MatElem) = nrows(M) == 0 ? 0 : rank(M)
form_rank(M::FloatGram) = size(M, 1) == 0 ? 0 : LinearAlgebra.rank(M)

_n_of(M::MatElem) = nrows(M)
_n_of(M::FloatGram) = size(M, 1)

# --- M2: countPosDiagEntries / countNegDiagEntries (internal; QQ and RR only) ------

function _pos_neg_diag_counts(M::Union{MatElem, Matrix{Float64}})
    M isa MatElem && !(base_ring(M) isa QQField) &&
        error("Only implemented over QQ and RR")
    _is_symmetric_square(M) || error("Matrix is not symmetric")
    A = _is_diag(M) ? M : diagonalize_via_congruence(M)
    pos = 0
    neg = 0
    for i in 1:_n_of(A)
        A[i, i] > 0 && (pos += 1)
        A[i, i] < 0 && (neg += 1)
    end
    (pos, neg)
end
_pos_neg_diag_counts(::Any) = error("Only implemented over QQ and RR")

# --- M2: getSignature (QQ and RR only) ----------------------------------------------
# Over QQ, reuses Hecke's signature_tuple on the quadratic space (verified identical
# to M2's diagonal-count signature on a 200-form sweep). The RR (float) path keeps
# M2's diagonal counts.

"""
    form_signature(beta::GWClass)

The signature of a symmetric bilinear form over ``\\mathbb{Q}`` or ``\\mathbb{R}``: after diagonalizing,
the number of positive diagonal entries minus the number of negative ones.
Together with the rank it classifies forms over ``\\mathbb{R}``; over ``\\mathbb{Q}`` it is one of the
invariants entering [`is_isomorphic_form`](@ref).

# Examples
```julia-repl
julia> form_signature(GWClass([0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]))
1

julia> form_signature(diagonal_form(QQ, (1, -1, 1)))
1
```

See also [`form_rank`](@ref), [`integral_discriminant`](@ref),
[`is_isomorphic_form`](@ref).
"""
function form_signature(beta::GWClass)
    M = gw_matrix(beta)
    if M isa MatElem && base_ring(M) isa QQField
        st = signature_tuple(quadratic_space(QQ, M))   # (pos, zero, neg)
        return st[1] - st[3]
    end
    pos, neg = _pos_neg_diag_counts(M)
    pos - neg
end

# --- M2: getIntegralDiscriminant (QQ only) ------------------------------------------

"""
    integral_discriminant(beta::GWClass)

A squarefree integral representative of the discriminant of a form over ``\\mathbb{Q}``:
the square class of the determinant of any Gram matrix representing `beta`,
normalized to a squarefree integer. The discriminant is one of the invariants
classifying rational forms (see [`is_isomorphic_form`](@ref)).

# Examples
```julia-repl
julia> beta = GWClass(QQ[1 4 7; 4 3 -1; 7 -1 5]);

julia> integral_discriminant(beta)
-269
```

See also [`form_signature`](@ref), [`hasse_witt_invariant`](@ref).
"""
function integral_discriminant(beta::GWClass)
    gw_base_field(beta) == QQ || error("GrothendieckWittClass is not over QQ")
    _squarefree_part(prod(diagonal_entries(beta); init = one(QQ)))
end

# --- M2: getRelevantPrimes (QQ only) ------------------------------------------------
# Primes dividing the entries of the squarefree diagonal representative, in
# first-seen order with each entry's factors ascending (M2's order — Test 24
# expects [23, 5, 47]).

"""
    relevant_primes(beta::GWClass)

A finite list of primes containing every prime at which the Hasse–Witt
invariant of the rational form `beta` can be nontrivial. The Hasse–Witt
invariants of a form equal 1 at all but finitely many primes ([S73, IV
§3.3]); since they are products of Hilbert symbols of the diagonal entries, it
suffices to take the primes dividing the entries of a squarefree diagonal
representative.

# Examples
```julia-repl
julia> relevant_primes(diagonal_form(QQ, (6, 7, 22)))
4-element Vector{ZZRingElem}:
 2
 3
 7
 11
```

# References
- [S73] J. P. Serre, *A course in arithmetic*, Springer-Verlag, 1973.

See also [`hasse_witt_invariant`](@ref).
"""
function relevant_primes(beta::GWClass)
    gw_base_field(beta) == QQ || error("GrothendieckWittClass is not over QQ")
    D = diagonal_entries(diagonal_class(beta))
    L = ZZRingElem[]
    for x in D
        for p in _prime_factors(ZZ(x))
            p in L || push!(L, p)
        end
    end
    L
end

# --- M2: getHasseWittInvariant ------------------------------------------------------
# Product of pairwise Hilbert symbols over i < j (M2's convention; some sources
# use i <= j — do not change this).

"""
    hasse_witt_invariant(beta::GWClass, p)
    hasse_witt_invariant(L::AbstractVector, p)

The Hasse–Witt invariant at the prime `p` of a form over ``\\mathbb{Q}``: for a
diagonalization ``⟨a_1, …, a_n⟩``, the product ``∏_{i<j} (a_i, a_j)_p`` of
pairwise Hilbert symbols (see [`hilbert_symbol_padic`](@ref)). The list
variant takes the diagonal entries directly.

The invariant equals 1 for all but finitely many primes — for `p` not
dividing any entry of a squarefree diagonal representative it is
automatically 1 — so only the [`relevant_primes`](@ref) need checking.

# Examples
```julia-repl
julia> beta = GWClass(QQ[1 4 7; 4 3 -1; 7 -1 5]);

julia> hasse_witt_invariant(beta, 7)
1

julia> hasse_witt_invariant([6, 7, 22], 2)
-1
```

See also [`hilbert_symbol_padic`](@ref), [`relevant_primes`](@ref),
[`is_isomorphic_form`](@ref).
"""
function hasse_witt_invariant(L::AbstractVector, p)
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    f = [_squarefree_part(x) for x in L]
    a = 1
    for i in 1:(length(f) - 1), j in (i + 1):length(f)
        a *= hilbert_symbol_padic(f[i], f[j], p)
    end
    a
end

# The class-level variant reuses Hecke's hasse_invariant (∏_{i<j} convention,
# verified identical to M2's on a 200-form sweep at 2/3/all relevant primes).
# The list variant above stays as M2's own product so that arbitrary
# user-supplied diagonal lists behave exactly as in M2.
function hasse_witt_invariant(beta::GWClass, p)
    gw_base_field(beta) == QQ ||
        error("method is only implemented over the rational numbers")
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    hasse_invariant(quadratic_space(QQ, gw_matrix(beta)), ZZ(p))
end
