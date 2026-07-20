# Port of M2 Code/MatrixMethods.m2.
#
# M2's dispatcher sends fields (and étale domains) to diagonalizeViaCongruenceField
# and other rings to the fraction-free diagonalizeViaCongruenceRing.
#
# The float (ℝ/ℂ) methods never share code paths that would push tolerances into
# the exact ℚ/𝔽_q arithmetic. Note that M2 itself uses exact `== 0` pivot tests
# over RR/CC in this routine — ported as-is.

# M2: isDiagonal / isSquareAndSymmetric — reuse Oscar's is_diagonal/is_symmetric
# for exact matrices and LinearAlgebra's isdiag/issymmetric for floats
# (issymmetric on Complex checks transpose-symmetry, not hermitian — correct here)
_is_diag(A::MatElem) = is_diagonal(A)
_is_diag(A::FloatGram) = LinearAlgebra.isdiag(A)

_is_symmetric_square(M::MatElem) = nrows(M) == ncols(M) && is_symmetric(M)
_is_symmetric_square(M::FloatGram) =
    size(M, 1) == size(M, 2) && LinearAlgebra.issymmetric(M)

# Symmetric Gaussian elimination by congruence (M2: diagonalizeViaCongruenceField).
# Mutates A; A must be square symmetric over a field. `divide(t, d)` and the
# index/arithmetic primitives are shared shapes, but each method keeps its own
# arithmetic (exact vs float).
function _diagonalize_field!(A, n, divide)
    _is_diag(A) && return A     # M2 short-circuits diagonal input
    for col in 1:n
        if iszero(A[col, col])
            for row in (col + 1):n
                if !iszero(A[row, col])
                    if iszero(A[row, row])
                        # row operation + matching column operation
                        for k in 1:n
                            A[col, k] += A[row, k]
                        end
                        for k in 1:n
                            A[k, col] += A[k, row]
                        end
                    else
                        # row and column swaps
                        for k in 1:n
                            A[col, k], A[row, k] = A[row, k], A[col, k]
                        end
                        for k in 1:n
                            A[k, col], A[k, row] = A[k, row], A[k, col]
                        end
                    end
                    break
                end
            end
        end
        if !iszero(A[col, col])
            for row in (col + 1):n
                c = -divide(A[row, col], A[col, col])
                for k in 1:n
                    A[row, k] += c * A[col, k]
                end
                for k in 1:n
                    A[k, row] += c * A[k, col]
                end
            end
        end
    end
    A
end

# M2: diagonalizeViaCongruence (dispatcher). Fields and étale DOMAINS (M2: zero-dim
# quotient with prime zero ideal, via toField) take the field path with division;
# any other ring takes the fraction-free ring path (diagonalizeViaCongruenceRing).
# Oscar's divexact on MPolyQuoRingElem provides the division over étale domains.
"""
    diagonalize_via_congruence(M)

A diagonal matrix congruent to the symmetric matrix `M`, computed by symmetric
Gaussian elimination (simultaneous row and column operations). `M` may be
defined over a field, a finite étale algebra, or as a `Matrix{Float64}` /
`Matrix{ComplexF64}` for ``\\mathbb{R}`` / ``\\mathbb{C}``. The order (and scaling) of the diagonal
entries is an artifact of the algorithm and is not normalized; use
[`diagonal_class`](@ref) for square-class-reduced entries.

Over étale algebras that are not domains the elimination is performed
fraction-free (rows and columns are scaled by the pivot instead of divided),
so the result is congruent but its entries may carry square factors.

# Examples
```julia-repl
julia> diagonalize_via_congruence(QQ[0 2; 2 0])
[4    0]
[0   -1]

julia> S, (y,) = polynomial_ring(QQ, ["y"]);

julia> A, _ = quo(S, ideal(S, [y^2 - 1]));

julia> diagonalize_via_congruence(A[A(1) A(2); A(2) A(y)])
[1       0]
[0   y - 4]
```

See also [`diagonal_class`](@ref), [`diagonal_entries`](@ref).
"""
function diagonalize_via_congruence(M::MatElem)
    _is_symmetric_square(M) || error("matrix is not symmetric")
    R = base_ring(M)
    R isa Field && return _diagonalize_field!(deepcopy(M), nrows(M), divexact)
    if R isa MPolyQuoRing && coefficient_ring(base_ring(R)) isa Field &&
       dim(modulus(R)) == 0 && is_prime(ideal(R, [zero(R)]))
        return _diagonalize_field!(deepcopy(M), nrows(M), divexact)
    end
    _diagonalize_ring!(deepcopy(M), nrows(M))
end

# M2: diagonalizeViaCongruenceRing — fraction-free symmetric elimination (no
# division: rows/columns are scaled by the pivot before clearing). No diagonal
# short-circuit, exactly as in M2.
function _diagonalize_ring!(A, n)
    for col in 1:n
        if iszero(A[col, col])
            for row in (col + 1):n
                if !iszero(A[row, col])
                    if iszero(A[row, row])
                        for k in 1:n
                            A[col, k] += A[row, k]
                        end
                        for k in 1:n
                            A[k, col] += A[k, row]
                        end
                    else
                        for k in 1:n
                            A[col, k], A[row, k] = A[row, k], A[col, k]
                        end
                        for k in 1:n
                            A[k, col], A[k, row] = A[k, row], A[k, col]
                        end
                    end
                    break
                end
            end
        end
        if !iszero(A[col, col])
            for row in (col + 1):n
                t = A[row, col]
                piv = A[col, col]
                for k in 1:n                       # M2: rowMult(A, row, piv)
                    A[row, k] *= piv
                end
                for k in 1:n                       # M2: columnMult(A, row, piv)
                    A[k, row] *= piv
                end
                for k in 1:n                       # M2: rowAdd(A, row, -t, col)
                    A[row, k] -= t * A[col, k]
                end
                for k in 1:n                       # M2: columnAdd(A, row, -t, col)
                    A[k, row] -= t * A[k, col]
                end
            end
        end
    end
    A
end

function diagonalize_via_congruence(M::Matrix{T}) where {T <: Union{Float64, ComplexF64}}
    _is_symmetric_square(M) || error("matrix is not symmetric")
    A = _diagonalize_field!(copy(M), size(M, 1), /)
    # M2 zeroes residual off-diagonal entries over inexact fields:
    # A = diagonalMatrix apply(n, i -> A_(i,i))
    LinearAlgebra.diagm(0 => LinearAlgebra.diag(A))
end

# --- per-field simplification (M2: diagonalizeAndSimplifyViaCongruence) ----------

function _diagonalize_and_simplify(M::MatElem)
    _is_symmetric_square(M) || error("matrix is not symmetric")
    A = diagonalize_via_congruence(M)
    R = base_ring(M)
    n = nrows(A)
    if R isa QQField
        # squarefree part of each diagonal entry
        for i in 1:n
            A[i, i] = R(_squarefree_part(A[i, i]))
        end
    elseif R isa FinField && characteristic(R) != 2
        # squares -> 1; nonsquares -> a fixed nonsquare representative, which is -1
        # unless -1 is itself a square, in which case it is the first nonsquare on
        # the diagonal
        non_square_rep = R(-1)
        if is_square(R(-1))
            for i in 1:n
                if !iszero(A[i, i]) && !is_square(A[i, i])
                    non_square_rep = A[i, i]
                    break
                end
            end
        end
        for i in 1:n
            if !iszero(A[i, i])
                A[i, i] = is_square(A[i, i]) ? one(R) : non_square_rep
            end
        end
    end
    # other fields: M2 (v2.0) applies no simplification and returns the
    # diagonalization unchanged
    A
end

function _diagonalize_and_simplify(M::Matrix{Float64})
    _is_symmetric_square(M) || error("matrix is not symmetric")
    A = diagonalize_via_congruence(M)
    n = size(A, 1)
    for i in 1:n
        # M2 over RR: positive entries -> 1, negative entries -> -1
        A[i, i] > 0 && (A[i, i] = 1.0)
        A[i, i] < 0 && (A[i, i] = -1.0)
    end
    A
end

function _diagonalize_and_simplify(M::Matrix{ComplexF64})
    _is_symmetric_square(M) || error("matrix is not symmetric")
    A = diagonalize_via_congruence(M)
    n = size(A, 1)
    for i in 1:n
        # M2 over CC: every nonzero entry -> 1
        iszero(A[i, i]) || (A[i, i] = 1.0 + 0.0im)
    end
    A
end

# --- nondegenerate part (M2: getNondegeneratePartDiagonal, internal) --------------

function _nondegenerate_part_diagonal(M::MatElem)
    A = diagonalize_via_congruence(M)
    keep = [i for i in 1:nrows(A) if !iszero(A[i, i])]
    diagonal_matrix(base_ring(M), [A[i, i] for i in keep])
end

function _nondegenerate_part_diagonal(M::Matrix{T}) where {T <: Union{Float64, ComplexF64}}
    A = diagonalize_via_congruence(M)
    keep = [i for i in 1:size(A, 1) if !iszero(A[i, i])]
    LinearAlgebra.diagm(0 => T[A[i, i] for i in keep])
end
