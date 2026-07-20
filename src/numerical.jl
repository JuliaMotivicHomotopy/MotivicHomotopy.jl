# The ℂ numerical backend, on HomotopyContinuation.jl. Numerical computation is
# needed ONLY for the unstable degrees over ℂ: there the class carries a scalar
# (a k^× factor), so the actual complex roots must be found. These add methods to
# the unstable-degree functions for floating-point input; HC names are qualified
# `HC.` so they never collide with Oscar's exports.
#
# The stable degrees do NOT use this backend. Over ℂ a symmetric bilinear form is
# determined by its rank, and the stable A¹-degree is the identity form of rank
# equal to the dimension of the (global or local) coordinate algebra — a discrete
# invariant obtained from the exact computation over ℚ, with no need to solve for
# roots. See `global_A1_degree` / `local_A1_degree` in degrees.jl.
#
# The unstable methods below find univariate roots via the companion matrix,
# cancel roots of numerator and denominator that agree within `linear_tolerance`,
# and form the scalar from the root-product resultant.

# ---------- univariate polynomial machinery ----------------------------------------

_to_c64(c) = c isa HC.Expression ?
    ComplexF64(HC.ModelKit.to_number(c)) : ComplexF64(c)

# ascending coefficient vector of a univariate expression
function _poly_coeffs(f::HC.Expression, var)
    E, C = HC.exponents_coefficients(f, [var])
    d = maximum(E)
    cv = zeros(ComplexF64, d + 1)
    for (k, e) in enumerate(E[1, :])
        cv[e + 1] = _to_c64(C[k])
    end
    cv
end

# degree-many roots, via the companion matrix
function _poly_roots(f::HC.Expression, var)
    cv = _poly_coeffs(f, var)
    d = length(cv) - 1
    d <= 0 && return ComplexF64[]
    C = zeros(ComplexF64, d, d)
    for i in 2:d
        C[i, i - 1] = 1
    end
    for i in 1:d
        C[i, d] = -cv[i] / cv[d + 1]
    end
    LinearAlgebra.eigvals(C)
end

_lead_coefficient(f::HC.Expression, var) = _poly_coeffs(f, var)[end]

function _the_variable(f::HC.Expression, g::HC.Expression)
    vars = unique!(vcat(HC.variables(f), HC.variables(g)))
    length(vars) == 1 ||
        error("both input polynomials must be defined over the same univariate polynomial ring")
    vars[1]
end

# cancel roots of L1 against L2 within tol, pairwise
function _remove_common_approx(L1, L2, tol)
    L1new = ComplexF64[]
    L2new = copy(L2)
    for r in L1
        i = findfirst(s -> abs(r - s) < tol, L2new)
        if i === nothing
            push!(L1new, r)
        else
            deleteat!(L2new, i)
        end
    end
    (L1new, L2new)
end

# ---------- unstable degrees over CC ------------------------------------------------

function global_unstable_A1_degree(f::HC.Expression, g::HC.Expression;
                                   linear_tolerance::Real = 1e-6)
    linear_tolerance > 0 || error("linearTolerance must be a positive number")
    var = _the_variable(f, g)
    lcg = _lead_coefficient(g, var)
    r1, r2 = _remove_common_approx(_poly_roots(f, var), _poly_roots(g, var),
                                   linear_tolerance)
    d = length(r1)
    d > length(r2) || error("the rational function is not pointed after reduction")
    # (-1)^((d^2-d)/2) * resultant(fr, gr), with fr monic from the roots
    scalar = ComplexF64((-1)^div(d^2 - d, 2))
    for a in r1
        scalar *= lcg * prod((a - b for b in r2); init = one(ComplexF64))
    end
    GWuClass(Matrix{ComplexF64}(LinearAlgebra.I, d, d), scalar)
end

function local_unstable_A1_degree(f::HC.Expression, g::HC.Expression, r::Number;
                                  linear_tolerance::Real = 1e-6)
    linear_tolerance > 0 || error("linearTolerance must be a positive number")
    var = _the_variable(f, g)
    lcg = _lead_coefficient(g, var)
    r1, r2 = _remove_common_approx(_poly_roots(f, var), _poly_roots(g, var),
                                   linear_tolerance)
    d = length(r1)
    d > length(r2) || error("the rational function is not pointed after reduction")
    rc = ComplexF64(r)
    frr = prod((rc - a for a in r1); init = one(ComplexF64))
    abs(frr) > linear_tolerance &&
        error("the field element is not a zero of the function after reduction")
    k = count(a -> abs(a - rc) < linear_tolerance, r1)
    roots_not_r = filter(a -> abs(a - rc) >= linear_tolerance, r1)
    ld_denom = prod((rc - a for a in roots_not_r); init = one(ComplexF64))
    grr = lcg * prod((rc - b for b in r2); init = one(ComplexF64))
    scalar = ComplexF64((-1)^div(k^2 - k, 2)) * (grr / ld_denom)^k
    GWuClass(Matrix{ComplexF64}(LinearAlgebra.I, k, k), scalar)
end
