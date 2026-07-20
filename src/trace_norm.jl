# Ports of getMultiplicationMatrix / getTrace / getNorm (M2: Code/TraceAndNorm.m2)
# and isFiniteEtaleAlgebra (M2: Code/GrothendieckWittClasses.m2).

# M2: getMultiplicationMatrix(C, a) — the matrix over the coefficient field of
# multiplication by a on a vector-space basis of the zero-dimensional algebra C.
"""
    multiplication_matrix(A, a)
    multiplication_matrix(S, I, b)

The matrix, over the coefficient field ``k``, of multiplication by an element
on a monomial basis of a finite-dimensional ``k``-algebra. The algebra is
given either directly as a quotient ring `A` (an `MPolyQuoRing`) with `a` an
element coercible into it, or as a polynomial ring `S` with an ideal `I` and
`b` an element of `S` (the algebra then being ``S/I``). The basis is the
standard monomials in ascending order.

# Examples
```julia-repl
julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(S, [x^2 + y^2 + 1, 3*x + 2]);

julia> A, _ = quo(S, I);

julia> multiplication_matrix(A, 1 + y*x^2)
[   1   -52//81]
[4//9         1]
```

See also [`algebra_trace`](@ref), [`algebra_norm`](@ref).
"""
function multiplication_matrix(A::MPolyQuoRing, a)
    aA = A(a)                      # M2: isRingElement — errors if not coercible
    # M2's basis lists standard monomials ascending; Oscar's monomial_basis is
    # descending — reversed so the matrix matches M2's literally (Tests 38/39)
    B = reverse(monomial_basis(A))
    n = length(B)
    kk = coefficient_ring(base_ring(A))
    M = zero_matrix(kk, n, n)
    for j in 1:n
        f = lift(simplify(aA * A(B[j])))   # reduced representative, in span of B
        for i in 1:n
            M[i, j] = coeff(f, B[i])
        end
    end
    M
end

# M2: getMultiplicationMatrix(S, I, b) — same, for the quotient S/I.
multiplication_matrix(S::MPolyRing, I::MPolyIdeal, b) =
    multiplication_matrix(quo(S, I)[1], b)

# M2: getTrace / getNorm
"""
    algebra_trace(A, a)
    algebra_trace(S, I, b)

The trace over ``k`` of an element of a finite-dimensional ``k``-algebra:
the trace of its [`multiplication_matrix`](@ref). Accepts the same input
shapes as `multiplication_matrix` (a quotient ring and an element, or a
polynomial ring, ideal, and element).

# Examples
```julia-repl
julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(S, [x^2 + y^2 + 1, 3*x + 2]);

julia> algebra_trace(S, I, 1 + y*x^2)
2
```

See also [`multiplication_matrix`](@ref), [`algebra_norm`](@ref),
[`transfer_gw`](@ref).
"""
algebra_trace(args...) = tr(multiplication_matrix(args...))

"""
    algebra_norm(A, a)
    algebra_norm(S, I, b)

The norm over ``k`` of an element of a finite-dimensional ``k``-algebra: the
determinant of its [`multiplication_matrix`](@ref). Accepts the same input
shapes as `multiplication_matrix`.

# Examples
```julia-repl
julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(S, [x^2 + y^2 + 1, 3*x + 2]);

julia> A, _ = quo(S, I);

julia> algebra_norm(A, 1 + y*x^2)
937//729
```

See also [`multiplication_matrix`](@ref), [`algebra_trace`](@ref).
"""
algebra_norm(args...) = det(multiplication_matrix(args...))

# M2: isFiniteEtaleAlgebra — a zero-dimensional algebra over a field whose trace
# form is nondegenerate. Fields short-circuit in M2 before this is consulted.
is_finite_etale_algebra(::Field) = true
function is_finite_etale_algebra(A::MPolyQuoRing)
    kk = coefficient_ring(base_ring(A))
    kk isa Field || return false
    dim(modulus(A)) == 0 || return false
    B = monomial_basis(A)
    n = length(B)
    T = zero_matrix(kk, n, n)
    for i in 1:n, j in 1:n
        T[i, j] = algebra_trace(A, A(B[i]) * A(B[j]))
    end
    !iszero(det(T))
end
is_finite_etale_algebra(::Any) = false
