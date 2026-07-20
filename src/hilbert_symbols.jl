# Port of M2 Code/HilbertSymbols.m2.
#
# The p-adic symbol reuses Hecke's `hilbert_symbol`, verified identical to M2's
# getHilbertSymbol (Serre III conventions, getSquareSymbol-based) on a 25,200-case
# sweep over p in {2,...,17}, a,b in ±[1,30], plus all M2 Test 25 values. M2's
# getSquareSymbol is therefore not ported. M2's input validation is reproduced
# here.

const _RatLike = Union{Integer, Rational, ZZRingElem, QQFieldElem}

# M2: getHilbertSymbol(a, b, p) — the Hilbert symbol (a,b)_p over QQ_p
"""
    hilbert_symbol_padic(a, b, p)

The Hilbert symbol ``(a, b)_p`` of two nonzero rational numbers, viewed as
elements of ``\\mathbb{Q}_p``:

``(a, b)_p = 1`` if ``z^2 = ax^2 + by^2`` has a nonzero solution over
``\\mathbb{Q}_p``, and ``-1`` otherwise ([S73, Chapter III]).

Products of Hilbert symbols compute the [`hasse_witt_invariant`](@ref), a key
step in classifying rational forms and certifying their (an)isotropy.

The name carries the `_padic` suffix because Oscar itself exports
`hilbert_symbol` (which this function calls internally).

# Examples

``z^2 = 2x^2 + y^2`` has the solution ``(1, 0, 3)`` mod 7, hence a 7-adic
solution by Hensel's lemma, while ``z^2 = 7x^2 + 3y^2`` has no nonzero
solution mod 7:

```julia-repl
julia> hilbert_symbol_padic(2, 1, 7)
1

julia> hilbert_symbol_padic(7, 3, 7)
-1

julia> hilbert_symbol_padic(2, 2, 2)
1

julia> hilbert_symbol_padic(2, 3, 2)
-1
```

# References
- [S73] J. P. Serre, *A course in arithmetic*, Springer-Verlag, 1973.

See also [`hilbert_symbol_real`](@ref), [`hasse_witt_invariant`](@ref).
"""
function hilbert_symbol_padic(a::_RatLike, b::_RatLike, p::Union{Integer, ZZRingElem})
    (iszero(a) || iszero(b)) &&
        error("first two arguments of getHilbertSymbol must be nonzero")
    is_prime(ZZ(p)) || error("third argument of getHilbertSymbol must be prime")
    hilbert_symbol(QQ(a), QQ(b), ZZ(p))
end

# M2: getHilbertSymbolReal(a, b) — the symbol (a,b)_RR: -1 iff both are negative
"""
    hilbert_symbol_real(a, b)

The Hilbert symbol ``(a, b)_{\\mathbb{R}}`` of two nonzero rational numbers viewed as
real numbers: ``-1`` if ``z^2 = ax^2 + by^2`` has no nonzero real solution —
which happens exactly when both `a` and `b` are negative — and ``1``
otherwise ([S73, Chapter III]).

# Examples
```julia-repl
julia> hilbert_symbol_real(-3, -2//3)
-1

julia> hilbert_symbol_real(3, -5)
1
```

# References
- [S73] J. P. Serre, *A course in arithmetic*, Springer-Verlag, 1973.

See also [`hilbert_symbol_padic`](@ref), [`form_signature`](@ref).
"""
function hilbert_symbol_real(a::_RatLike, b::_RatLike)
    (iszero(a) || iszero(b)) &&
        error("the arguments of getHilbertSymbolReal must be nonzero")
    (a < 0 && b < 0) ? -1 : 1
end
