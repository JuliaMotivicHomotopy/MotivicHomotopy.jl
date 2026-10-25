# Port of M2 Code/GWTransfer.m2: transferGW — the transfer GW(L) -> GW(k) along a
# finite étale algebra L/k, induced by the trace Tr_{L/k}: L -> k. The transferred
# Gram matrix is built from the Gram matrix over L itself, not from a
# diagonalization, so nothing has to be inverted in L, which need not be a field.
# (The GWClass constructor has already checked that L is finite étale.)

"""
    transfer_gw(beta::GWClass)

The image of a Grothendieck–Witt class over a finite étale algebra ``L/k``
under the transfer map ``\\text{GW}(L) → \\text{GW}(k)`` induced by the trace
``\\text{Tr}_{L/k} : L → k``. A form ``β`` on ``L^n`` is sent to the form
``(v, w) \\mapsto \\text{Tr}_{L/k}(β(v, w))`` on ``L^n`` regarded as a
``k``-vector space, so the rank is multiplied by ``[L : k]``. For example,
``⟨a⟩`` is sent to the scaled trace form ``(x, y) \\mapsto \\text{Tr}_{L/k}(axy)``,
of rank ``[L : k]``.

The Gram matrix is written in the ``k``-basis ``e_i b_s`` of ``L^n``, with ``b``
running fastest, where ``e_1, …, e_n`` is the standard basis of ``L^n`` and
``b_1, …, b_d`` are the standard monomials of ``L`` in the order used by
[`multiplication_matrix`](@ref): its ``(i, j)`` block is
``(\\text{Tr}_{L/k}(β_{ij} b_s b_t))_{s,t}``. Its determinant is
``\\det(T)^n N_{L/k}(\\det β)``, where ``T = (\\text{Tr}_{L/k}(b_s b_t))_{s,t}`` is
the trace form of ``L/k``.

!!! note
    The transferred form is nondegenerate exactly when ``\\det β`` is a unit of
    ``L``, which always holds when ``L`` is a field. When ``L`` is a product of
    fields and ``\\det β`` is a zero divisor, ``β`` is degenerate on a factor of
    ``L``, and `transfer_gw` errors.

# Examples

Over ``L = \\mathbb{Q}[t]/(t^2 - 1) \\cong \\mathbb{Q} \\times \\mathbb{Q}`` the
transfer is the sum of the specializations at ``t = 1`` and ``t = -1``:

```julia-repl
julia> S, (t,) = polynomial_ring(QQ, ["t"]);

julia> A, _ = quo(S, ideal(S, [t^2 - 1]));

julia> beta = GWClass(A[A(1) A(2); A(2) A(t)]);

julia> transfer_gw(beta)
[2   0   4   0]
[0   2   0   4]
[4   0   0   2]
[0   4   2   0]

julia> is_isomorphic_form(transfer_gw(beta), gw_direct_sum(GWClass(QQ[1 2; 2 1]), GWClass(QQ[1 2; 2 -1])))
true
```

See also [`algebra_trace`](@ref), [`algebra_norm`](@ref), [`GWClass`](@ref).
"""
function transfer_gw(alpha::GWClass)
    Alg = gw_algebra(alpha)
    Alg isa MPolyQuoRing ||
        error("transferGW requires a class over a finite étale algebra")
    kk = coefficient_ring(base_ring(Alg))
    M = gw_matrix(alpha)
    n = nrows(M)
    # det M is a unit exactly when multiplication by it is invertible
    iszero(algebra_norm(Alg, det(M))) &&
        error("the determinant of the form is a zero divisor in the algebra, " *
              "so its transfer is degenerate")
    mons = reverse(monomial_basis(Alg))   # the basis of multiplication_matrix
    d = length(mons)
    T = matrix(kk, [algebra_trace(Alg, mons[s] * mons[t]) for s in 1:d, t in 1:d])
    G = zero_matrix(kk, n * d, n * d)
    for j in 1:n, i in 1:j
        # Tr(M[i, j] b_s b_t) is entry (s, t) of T times multiplication by M[i, j]
        TM = T * multiplication_matrix(Alg, M[i, j])
        for s in 1:d, t in 1:d
            G[(i - 1) * d + s, (j - 1) * d + t] = TM[s, t]
            G[(j - 1) * d + t, (i - 1) * d + s] = TM[s, t]
        end
    end
    GWClass(G)
end
