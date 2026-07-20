# Port of M2 Code/Isotropy.m2: isAnisotropic / isIsotropic.
# Renamed is_anisotropic_form / is_isotropic_form: Oscar exports is_isotropic, and
# the pair is kept symmetric.

"""
    is_anisotropic_form(beta::GWClass)
    is_anisotropic_form(A)

Whether a symmetric bilinear form over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite field of odd
characteristic is anisotropic, i.e. has no nonzero vector ``v`` with
``β(v, v) = 0``. Computed as the statement that the anisotropic dimension
equals the dimension of the form.

What this takes per field: over ``\\mathbb{C}`` only rank-one forms are anisotropic; over
``\\mathbb{R}`` a form is anisotropic iff its diagonal entries are all positive or all
negative; over ``\\mathbb{Q}`` the Hasse–Minkowski principle ([L05, VI.3.1]) reduces the
question to the completions (forms of rank ≥ 5 over ``\\mathbb{Q}_p`` are always
isotropic ([S73, IV Theorem 6]), so only finitely many invariant computations
are needed); over a finite field a nondegenerate form is anisotropic iff its
rank is ≤ 2 and it is not hyperbolic.

# Examples
```julia-repl
julia> is_anisotropic_form(GWClass(ComplexF64[2 0; 0 5]))
false

julia> is_anisotropic_form(GWClass([3.0 0 0; 0 5 0; 0 0 7]))
true

julia> is_anisotropic_form(GWClass(GF(7)[1 0 0; 0 1 0; 0 0 1]))
false
```

# References
- [L05] T. Y. Lam, *Introduction to quadratic forms over fields*, American Mathematical Society, 2005.
- [S73] J. P. Serre, *A course in arithmetic*, Springer-Verlag, 1973.

See also [`is_isotropic_form`](@ref), [`anisotropic_dimension`](@ref),
[`anisotropic_part`](@ref).
"""
function is_anisotropic_form(A::Union{MatElem, Matrix{Float64}, Matrix{ComplexF64}})
    if A isa MatElem
        R = base_ring(A)
        (R isa QQField || (R isa FinField && characteristic(R) != 2)) ||
            error(_UNSUPPORTED_FIELD_ERR)
    end
    _n_of(A) == anisotropic_dimension(A)
end
is_anisotropic_form(alpha::GWClass) = is_anisotropic_form(gw_matrix(alpha))

"""
    is_isotropic_form(beta::GWClass)
    is_isotropic_form(A)

Whether a symmetric bilinear form over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite field of odd
characteristic is isotropic — the negation of
[`is_anisotropic_form`](@ref); see there for the per-field criteria.

# Examples
```julia-repl
julia> is_isotropic_form(diagonal_form(QQ, (1, -1)))
true

julia> is_isotropic_form(GWClass(GF(7)[3 0; 0 3]))
false
```

See also [`gw_witt_index`](@ref), [`anisotropic_dimension`](@ref).
"""
is_isotropic_form(A::Union{MatElem, Matrix{Float64}, Matrix{ComplexF64}}) =
    !is_anisotropic_form(A)
is_isotropic_form(alpha::GWClass) = !is_anisotropic_form(alpha)
