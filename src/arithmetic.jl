# Ports of the arithmetic helpers from M2 Code/ArithmeticMethods.m2.
#
# M2 isGFSquare is not reimplemented: Oscar's `is_square` on finite-field elements
# is an exact equivalent.

# M2: getSquarefreePart — the smallest-magnitude integer in the square class of n.
function _squarefree_part(n::ZZRingElem)
    iszero(n) && return zero(ZZ)
    f = factor(abs(n))
    sign(n) * prod(p^(e % 2) for (p, e) in f; init = one(ZZ))
end
_squarefree_part(q::QQFieldElem) = _squarefree_part(numerator(q) * denominator(q))
_squarefree_part(n::Integer) = _squarefree_part(ZZ(n))
_squarefree_part(q::Rational) = _squarefree_part(QQ(q))

# M2: getPrimeFactors — sorted prime divisors of |n|; reuses Oscar's prime_divisors
# (sorted explicitly: M2 sorts, prime_divisors does not guarantee order)
_prime_factors(n::ZZRingElem) = sort!(prime_divisors(abs(n)))

# M2: getPadicValuation — reuses Oscar's `valuation` (verified equivalent on probe);
# M2's zero check is reproduced.
"""
    padic_valuation(a, p)

The ``p``-adic valuation of a nonzero integer or rational number `a`: the
integer ``n`` with ``a = u·p^n`` for a unit ``u`` of ``\\mathbb{Z}_p``. Errors on
``a = 0``.

# Examples

``363/7 = 3·11^2/7``, so the 11-adic valuation is 2:

```julia-repl
julia> padic_valuation(363//7, 11)
2
```
"""
function padic_valuation(n::Union{Integer, ZZRingElem}, p::Union{Integer, ZZRingElem})
    iszero(n) && error("Trying to find prime factorization of 0")
    valuation(ZZ(n), ZZ(p))
end
function padic_valuation(q::Union{Rational, QQFieldElem}, p::Union{Integer, ZZRingElem})
    iszero(q) && error("Trying to find prime factorization of 0")
    valuation(QQ(q), ZZ(p))
end

# M2: isPadicSquare — reuses Hecke's `is_local_square`, verified identical to M2's
# squarefree-part/valuation-parity/mod-8 logic on a 600-case sweep.
_is_padic_square(a, p) = is_local_square(QQ(a), ZZ(p))
