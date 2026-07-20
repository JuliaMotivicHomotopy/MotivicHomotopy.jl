# Port of M2 Code/SimplifiedRepresentatives.m2 (getDiagonalClass, getDiagonalEntries)
# plus the getDiagonalClass method for unstable classes from
# Code/UnstableGrothendieckWittClasses.m2.
#
# Cache: M2 stores the result under beta.cache.getDiagonalClass; we use
# beta.cache[:diagonal_class]. (M2's getSumDecomposition later OVERWRITES this
# same key — see decomposition.jl.) One deliberate difference: M2's first call
# returns a fresh twin of the cached object; we return the cached object itself
# — indistinguishable under ==.

# M2: getDiagonalClass GrothendieckWittClass
"""
    diagonal_class(beta)

A class isomorphic to `beta` with a diagonal Gram matrix, with simplified
diagonal entries: over ``\\mathbb{Q}`` each entry is replaced by its squarefree integral
representative; over a finite field entries become 1 or a fixed nonsquare;
over ``\\mathbb{R}`` entries become ±1 and over ``\\mathbb{C}`` they become 1. For a [`GWuClass`](@ref)
the stable part is diagonalized and the scalar carried through unchanged.

The result is cached on `beta`, so repeated calls are free. Note that
[`sum_decomposition`](@ref) *overwrites* this cache with its own
representative, so the value returned by `diagonal_class` can change after a
call to `sum_decomposition`.

# Examples
```julia-repl
julia> beta = GWClass(QQ[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5]);

julia> diagonal_class(beta)
[1    0     0     0]
[0   89     0     0]
[0    0   445     0]
[0    0     0   -55]
```

See also [`diagonalize_via_congruence`](@ref), [`diagonal_entries`](@ref).
"""
function diagonal_class(beta::GWClass)
    haskey(beta.cache, :diagonal_class) && return beta.cache[:diagonal_class]::GWClass
    result = GWClass(_diagonalize_and_simplify(gw_matrix(beta)))
    beta.cache[:diagonal_class] = result
    result
end

# M2: getDiagonalClass UnstableGrothendieckWittClass — simplifies the stable part,
# carries the scalar through unchanged.
function diagonal_class(beta::GWuClass)
    haskey(beta.cache, :diagonal_class) && return beta.cache[:diagonal_class]::GWuClass
    result = GWuClass(_diagonalize_and_simplify(gw_matrix(beta)), gw_scalar(beta))
    beta.cache[:diagonal_class] = result
    result
end

# M2: getDiagonalEntries — diagonalizes (without simplification) and reads off the
# diagonal. M2 installs this for GrothendieckWittClass only.
"""
    diagonal_entries(beta::GWClass)

The entries ``a_1, …, a_n`` such that ``β ≅ ⟨a_1, …, a_n⟩``, obtained by
diagonalizing the Gram matrix via congruence (without square-class
simplification) and reading off the diagonal. If `beta` is already diagonal
the entries are returned as-is.

# Examples
```julia-repl
julia> diagonal_entries(GWClass(QQ[3 0 0; 0 2 0; 0 0 7]))
3-element Vector{QQFieldElem}:
 3
 2
 7

julia> diagonal_entries(GWClass([0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]))
3-element Vector{Float64}:
  2.0
  1.0
 -0.5
```

See also [`diagonal_class`](@ref), [`diagonalize_via_congruence`](@ref).
"""
function diagonal_entries(beta::GWClass)
    M = diagonalize_via_congruence(gw_matrix(beta))
    _diag_of(M)
end

_diag_of(M::MatElem) = [M[i, i] for i in 1:nrows(M)]
_diag_of(M::Matrix{T}) where {T} = T[M[i, i] for i in 1:size(M, 1)]
