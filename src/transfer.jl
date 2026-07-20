# Port of M2 Code/GWTransfer.m2: transferGW — the map GW(L) -> GW(k) for a finite
# étale algebra L/k, computed by diagonalizing over L and applying the trace to
# each diagonal entry. (M2 performs no validation of the input algebra; the
# étale-domain check happens implicitly inside the diagonalization dispatch.)

"""
    transfer_gw(beta::GWClass)

The image of a Grothendieck–Witt class over a finite étale algebra ``L/k``
under the canonical transfer map ``\\text{GW}(L) → \\text{GW}(k)``, computed by
diagonalizing over ``L`` and applying the trace form: the result is the
diagonal form over ``k`` whose entries are the traces
([`algebra_trace`](@ref)) of the diagonal entries.

!!! note
    If the trace of a diagonal entry vanishes, the would-be output is
    degenerate and the constructor errors.

# Examples
```julia-repl
julia> S, (t,) = polynomial_ring(QQ, ["t"]);

julia> A, _ = quo(S, ideal(S, [t^2 - 1]));

julia> beta = GWClass(A[A(1) A(2); A(2) A(t)]);

julia> transfer_gw(beta)
[2    0]
[0   -8]
```

See also [`algebra_trace`](@ref), [`GWClass`](@ref).
"""
function transfer_gw(alpha::GWClass)
    Alg = gw_algebra(alpha)
    Alg isa MPolyQuoRing ||
        error("transferGW requires a class over a finite étale algebra")
    kk = coefficient_ring(base_ring(Alg))
    diagonal_form(kk, Tuple(algebra_trace(Alg, e) for e in diagonal_entries(alpha)))
end
