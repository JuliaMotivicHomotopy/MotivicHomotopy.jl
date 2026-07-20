# Port of M2 Code/Decomposition.m2: the QQ anisotropic-part reducers
# (Koprowski/Rothkegel Algorithms 5, 7, 8 over QQ), getAnisotropicPart,
# getSumDecomposition(Verbose), getSumDecompositionString.
#
# M2's solveCongruenceList/solveCongruencePair/computeExtendedEuclidean are not
# ported: Oscar's `crt` is the equivalent (representative may differ from M2's;
# downstream results agree as classes).

# M2: reduceAnisotropicPartQQDimension4 — <1> or <-1> by signature sign
function _reduce_aniso_dim4(beta::GWClass)
    anisotropic_dimension(beta) >= 4 ||
        error("anisotropic dimension of form is not >= 4")
    form_signature(beta) >= 0 ? diagonal_form(QQ, 1) : diagonal_form(QQ, -1)
end

# M2: reduceAnisotropicPartQQDimension3 — a <a> with aniso(beta + <-a>) == 2
function _reduce_aniso_dim3(beta::GWClass)
    anisotropic_dimension(beta) == 3 ||
        error("anisotropic dimension of form is not 3")
    d = integral_discriminant(beta)
    residues = ZZRingElem[]
    moduli = ZZRingElem[]
    for p in relevant_primes(beta)
        if isodd(padic_valuation(d, p))
            push!(residues, ZZ(d - 1))   # a unit ≡ d-1 mod p
            push!(moduli, p)
        else
            push!(residues, p)           # valuation exactly 1 at p
            push!(moduli, p^2)
        end
    end
    alpha = isempty(moduli) ? one(ZZ) : crt(residues, moduli)
    bigM = isempty(moduli) ? one(ZZ) : prod(moduli)
    # The congruences fix alpha only modulo bigM (and hence fix its square class
    # at every relevant finite prime); they say nothing about the real place, so
    # representatives alpha + k*bigM of opposite sign are all admissible at the
    # finite primes. Choose the sign for which beta ⊕ ⟨-alpha⟩ actually has
    # anisotropic dimension 2 — the split is valid only then. (A signature-based
    # heuristic is not enough: for signature ±1 either sign can over-reduce, and
    # which one does is representative-dependent.)
    r = mod(alpha, bigM)
    candidates = iszero(r) ? (bigM, -bigM) : (r, r - bigM)
    for cand in candidates
        a = _squarefree_part(cand)
        if anisotropic_dimension(gw_direct_sum(beta, diagonal_form(QQ, -a))) == 2
            return diagonal_form(QQ, a)
        end
    end
    diagonal_form(QQ, _squarefree_part(alpha))   # fallback (not expected)
end

# M2: getAnisotropicPartQQDimension2 — the anisotropic part when aniso dim == 2
function _aniso_part_qq_dim2(beta::GWClass)
    anisotropic_dimension(beta) == 2 ||
        error("anisotropic dimension of form is not 2")
    n = form_rank(beta)
    n == 2 && return beta                       # M2 shortcut

    # Step 1: make the Witt index 0 mod 4
    w = gw_witt_index(beta)
    q = beta
    if w % 4 != 0
        w = w % 4
        q = gw_direct_sum(q, hyperbolic_form(QQ, 2 * (4 - w)))
        n = n + 2 * (4 - w)
    end

    # Step 2: signed discriminant (Koprowski/Rothkegel convention)
    d = (-1)^div(n * (n - 1), 2) * integral_discriminant(q)

    # Step 3: relevant primes of the ORIGINAL beta, plus 2 (M2 appends 2 at the end)
    S = relevant_primes(beta)
    ZZ(2) in S || push!(S, ZZ(2))

    p = ZZ(2)
    local basisES, X
    while true
        s = length(S)
        # Step 5a: basis of the S-singular group: S ∪ {-1}
        basisES = vcat(S, ZZRingElem[ZZ(-1)])
        m = length(basisES)
        # Step 5c: Hasse-invariant exponent vector
        Wrows = Int[div(1 - hasse_witt_invariant(q, S[i]), 2) for i in 1:s]
        # Step 5b/5f: real place row when the discriminant is negative
        if d < 0
            abs(form_signature(q)) == 2 || error("getSignature isn't pm 2")
            Wrows = vcat(Int[form_signature(q) == 2 ? 0 : 1], Wrows)
        end
        # Step 5e: matrix of Hilbert-symbol exponents
        B = Int[div(1 - hilbert_symbol_padic(basisES[j], d, S[i]), 2)
                for i in 1:s, j in 1:m]
        # Step 5d: sign row when the discriminant is negative
        if d < 0
            B = vcat(reshape(Int[basisES[j] > 0 ? 0 : 1 for j in 1:m], 1, m), B)
        end
        # Step 5f: solve over GF(2)
        k2 = GF(2)
        ok, sol = can_solve_with_solution(matrix(k2, B),
                                          matrix(k2, length(Wrows), 1, Wrows);
                                          side = :right)
        if ok
            X = sol
            break
        end
        # enlarge S by the next prime not yet present (M2: nextPrime(p+1) loop)
        p = next_prime(p)
        while p in S
            p = next_prime(p)
        end
        push!(S, p)
    end

    alpha = one(ZZ)
    for j in eachindex(basisES)
        alpha *= basisES[j]^Int(lift(ZZ, X[j, 1]))
    end
    diagonal_form(QQ, (alpha, -_squarefree_part(alpha * d)))
end

# M2: getAnisotropicPartQQ — strip <a>'s until the dim-2/1 base cases
function _anisotropic_part_qq(beta::GWClass)
    beta = diagonal_class(beta)
    n = form_rank(beta)
    anisotropic_dimension(beta) == n && return beta

    output = diagonal_form(QQ, ())
    while anisotropic_dimension(beta) >= 4
        red = _reduce_aniso_dim4(beta)
        output = gw_direct_sum(output, red)
        a = gw_matrix(red)[1, 1]
        beta = gw_direct_sum(beta, diagonal_form(QQ, -a))
    end
    if anisotropic_dimension(beta) == 3
        red = _reduce_aniso_dim3(beta)
        output = gw_direct_sum(output, red)
        a = gw_matrix(red)[1, 1]
        beta = gw_direct_sum(beta, diagonal_form(QQ, -a))
    end
    if anisotropic_dimension(beta) == 2
        output = gw_direct_sum(output, _aniso_part_qq_dim2(beta))
    end
    if anisotropic_dimension(beta) == 1
        # n is the ORIGINAL rank: this branch is only reachable with beta untouched
        output = gw_direct_sum(output,
            diagonal_form(QQ, (-1)^div(n - 1, 2) * integral_discriminant(beta)))
    end
    output
end

# --- M2: getAnisotropicPart, per field ---------------------------------------------

"""
    anisotropic_part(beta::GWClass)
    anisotropic_part(A)

The anisotropic part of a symmetric bilinear form over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite
field of odd characteristic: the (unique up to isomorphism) anisotropic form
``β_a`` in the Witt decomposition ``β ≅ β_a ⊕ n·\\mathbb{H}``.

Over ``\\mathbb{C}``, ``\\mathbb{R}``, and finite fields this is a short computation from the rank,
signature, or discriminant. Over ``\\mathbb{Q}`` it uses the number-field algorithms of
Koprowski–Rothkegel [KR23]: ranks ≥ 4 are peeled off by signs of the
signature, rank 3 via a CRT-constructed splitting element, and the rank-2
base case via a Hilbert-symbol exponent system solved over GF(2).

# Examples
```julia-repl
julia> anisotropic_part(diagonal_form(QQ, (3, -3, 2, 5, 1, -9)))
[2   0]
[0   5]
```

# References
- [KR23] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic form over a number field*, Journal of Symbolic Computation, 2023.

See also [`anisotropic_dimension`](@ref), [`gw_witt_index`](@ref),
[`sum_decomposition`](@ref).
"""
function anisotropic_part(A::MatElem)
    R = base_ring(A)
    (R isa QQField || (R isa FinField && characteristic(R) != 2)) ||
        error(_UNSUPPORTED_FIELD_ERR)
    _is_symmetric_square(A) || error("Underlying matrix is not symmetric")
    if R isa QQField
        return gw_matrix(_anisotropic_part_qq(GWClass(_nondegenerate_part_diagonal(A))))
    end
    # finite field: <1> or <e> for dim 1; <1,-e> shape for dim 2
    diagA = diagonalize_via_congruence(A)
    r = form_rank(diagA)
    ad = anisotropic_dimension(A)
    if ad == 0
        return diagonal_matrix(R, elem_type(R)[])
    elseif ad == 1
        return matrix(R, 1, 1,
            [R((-1)^div(r - 1, 2)) * det(_nondegenerate_part_diagonal(diagA))])
    else
        return matrix(R, 2, 2,
            [one(R), zero(R), zero(R),
             R((-1)^div(r - 2, 2)) * det(_nondegenerate_part_diagonal(diagA))])
    end
end

function anisotropic_part(A::Matrix{Float64})
    _is_symmetric_square(A) || error("Underlying matrix is not symmetric")
    pos, neg = _pos_neg_diag_counts(diagonalize_via_congruence(A))
    pos > neg && return Matrix{Float64}(LinearAlgebra.I, pos - neg, pos - neg)
    pos < neg && return -Matrix{Float64}(LinearAlgebra.I, neg - pos, neg - pos)
    zeros(Float64, 0, 0)
end

function anisotropic_part(A::Matrix{ComplexF64})
    _is_symmetric_square(A) || error("Underlying matrix is not symmetric")
    anisotropic_dimension(A) == 0 ? zeros(ComplexF64, 0, 0) :
        ComplexF64[1.0 + 0.0im;;]
end

anisotropic_part(alpha::GWClass) = GWClass(anisotropic_part(gw_matrix(alpha)))

# --- M2 element printing for decomposition strings ---------------------------------
# M2 prints fp elements with balanced representatives (8 mod 13 -> "-5"), rationals
# as "a/b", and inexact ±1.0 as integers. Best-effort `string` elsewhere.

_m2_string(x::QQFieldElem) = isone(denominator(x)) ? string(numerator(x)) :
    string(numerator(x)) * "/" * string(denominator(x))
_m2_string(x::ZZRingElem) = string(x)
function _m2_string(x::FinFieldElem)
    v = lift(ZZ, x)
    p = order(parent(x))
    string(2 * v > p ? v - p : v)
end
_m2_string(x::Float64) = isinteger(x) ? string(Int(x)) : string(x)
_m2_string(x::ComplexF64) = (iszero(imag(x)) && isinteger(real(x))) ?
    string(Int(real(x))) : string(x)
_m2_string(x) = string(x)

# --- M2: getSumDecompositionVerbose -------------------------------------------------

function _sum_decomposition_verbose(beta::GWClass)
    kk = gw_base_field(beta)
    form_rank(beta) == 0 && return (GWClass(_diag_gram(kk, [])), "empty form")

    str = ""
    w = gw_witt_index(beta)
    w == 1 && (str *= "H")
    w > 1 && (str *= string(w) * "H")

    hyperbolic_part = hyperbolic_form(kk, 2 * w)
    alpha = anisotropic_part(beta)

    if form_rank(alpha) > 0
        for d in diagonal_entries(alpha)
            piece = "<" * _m2_string(d) * ">"
            str = isempty(str) ? piece : str * " + " * piece
        end
    end
    (gw_direct_sum(alpha, hyperbolic_part), str)
end

# --- M2: getSumDecomposition / getSumDecompositionString ---------------------------
# Cache: mirrors M2 exactly — the result OVERWRITES cache[:diagonal_class], so
# diagonal_class returns the sum decomposition afterwards.

"""
    sum_decomposition(beta)

A simplified diagonal representative of a [`GWClass`](@ref) or
[`GWuClass`](@ref) over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite field of odd characteristic:
the class rewritten as its [`anisotropic_part`](@ref) plus
[`gw_witt_index`](@ref)-many hyperbolic forms. For an unstable class the
decomposition is applied to the stable part and the scalar kept. Over ``\\mathbb{R}`` this
reflects the classification of a form by its rank and signature ([L05, II
Proposition 3.5]).

The result overwrites the `diagonal_class` cache slot on `beta`, so a later
[`diagonal_class`](@ref) call returns this representative.

# Examples
```julia-repl
julia> gamma = GWClass(QQ[1 2 3; 2 4 5; 3 5 6]);

julia> sum_decomposition(gamma)
[1   0    0]
[0   1    0]
[0   0   -1]

julia> delta = GWClass(GF(13)[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5]);

julia> sum_decomposition(delta)
[1   0   0    0]
[0   8   0    0]
[0   0   1    0]
[0   0   0   12]
```

# References
- [L05] T. Y. Lam, *Introduction to quadratic forms over fields*, American Mathematical Society, 2005.

See also [`sum_decomposition_string`](@ref), [`anisotropic_part`](@ref),
[`gw_witt_index`](@ref).
"""
function sum_decomposition(beta::GWClass)
    result = _sum_decomposition_verbose(beta)[1]
    beta.cache[:diagonal_class] = result
    result
end

"""
    sum_decomposition_string(beta)

A human-readable string for the [`sum_decomposition`](@ref) of a
[`GWClass`](@ref) or [`GWuClass`](@ref): hyperbolic summands are written `H`
(with a multiplicity prefix) and rank-one summands `<a>`. For an unstable
class the result is the pair `"(decomposition, scalar)"`.

# Examples
```julia-repl
julia> sum_decomposition_string(GWClass(QQ[1 2 3; 2 4 5; 3 5 6]))
"H + <1>"

julia> sum_decomposition_string(GWClass(GF(13)[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5]))
"H + <1> + <-5>"
```

See also [`sum_decomposition`](@ref).
"""
sum_decomposition_string(beta::GWClass) = _sum_decomposition_verbose(beta)[2]

function sum_decomposition(beta::GWuClass)
    simplified = sum_decomposition(stable_part(beta))
    result = GWuClass(simplified, gw_scalar(beta))
    beta.cache[:diagonal_class] = result
    result
end

sum_decomposition_string(beta::GWuClass) =
    "(" * sum_decomposition_string(stable_part(beta)) * ", " *
    _m2_string(gw_scalar(beta)) * ")"
