# Port of M2 Code/BuildingForms.m2: isCompatibleElement and the form builders
# makeDiagonalForm / makeHyperbolicForm / makePfisterForm and their unstable
# variants. (makeAntidiagonalUnstableForm lives with the unstable degrees,
# unstable_degrees.jl.)

# --- isCompatibleElement, as a ring-level coercion ------------------------------
# Returns the element coerced into kk, or `nothing` if incompatible.
# Float64 / ComplexF64 stand in for M2's RR / CC (InexactFieldFamily variants).

function _coerce_into(R::QQField, b)
    b isa Union{Integer, Rational, ZZRingElem, QQFieldElem} || return nothing
    R(b)
end

function _coerce_into(R::FinField, b)
    if b isa Union{Integer, ZZRingElem}
        return R(b)
    elseif b isa Union{Rational, QQFieldElem}
        q = QQ(b)
        den = R(denominator(q))
        iszero(den) && return nothing
        return R(numerator(q)) * inv(den)
    elseif b isa FinFieldElem
        parent(b) === R && return b
        (parent(b) isa FinField && order(parent(b)) == order(R)) || return nothing
        return try R(lift(ZZ, b)) catch; nothing end
    end
    nothing
end

function _coerce_into(R::MPolyQuoRing, b)
    (b isa MPolyQuoRingElem && parent(b) === R) || return nothing
    b
end

function _coerce_into(R::Field, b)   # other exact fields: element of R required
    (b isa RingElem && parent(b) === R) || return nothing
    b
end

_coerce_into(::Type{Float64}, b) = b isa Real ? Float64(b) : nothing
_coerce_into(::Type{ComplexF64}, b) = b isa Number ? ComplexF64(b) : nothing
_coerce_into(::Any, ::Any) = nothing

function _coerce_into_or_error(kk, b)
    c = _coerce_into(kk, b)
    c === nothing && error("scalar not compatible with field")   # M2's message
    c
end

# --- diagonal Gram construction per backend --------------------------------------

_diag_gram(kk::Ring, vals) = diagonal_matrix(kk, elem_type(kk)[v for v in vals])
_diag_gram(::Type{Float64}, vals) = LinearAlgebra.diagm(0 => Float64[v for v in vals])
_diag_gram(::Type{ComplexF64}, vals) = LinearAlgebra.diagm(0 => ComplexF64[v for v in vals])

# --- M2: makeDiagonalForm ---------------------------------------------------------

_diagonal_form(kk, entries) =
    GWClass(_diag_gram(kk, [_coerce_into_or_error(kk, b) for b in entries]))

"""
    diagonal_form(kk, (a₁, …, aₙ))
    diagonal_form(kk, a₁, a₂, …)

The [`GWClass`](@ref) of the diagonal form ``⟨a_1, …, a_n⟩`` over the field
or finite étale algebra `kk`: the block sum of the rank-one forms
``⟨a_i⟩ : k × k → k``, ``(x, y) ↦ a_i x y``. A single entry produces a
rank-one form. For forms over ``\\mathbb{R}`` / ``\\mathbb{C}`` pass `Float64` / `ComplexF64` as `kk`.

# Examples
```julia-repl
julia> diagonal_form(QQ, (3, 5, 7))
[3   0   0]
[0   5   0]
[0   0   7]

julia> diagonal_form(GF(29), 5//13)
[16]

julia> diagonal_form(Float64, 2)
1×1 Matrix{Float64}:
 2.0
```

See also [`hyperbolic_form`](@ref), [`pfister_form`](@ref),
[`diagonal_unstable_form`](@ref), [`diagonal_entries`](@ref).
"""
diagonal_form(kk, L::Tuple) = _diagonal_form(kk, collect(L))
diagonal_form(kk, entries...) = _diagonal_form(kk, collect(entries))

# --- M2: makeHyperbolicForm -------------------------------------------------------
# H = <1,-1>; rank n gives n/2 block copies, i.e. diag(1,-1,1,-1,...).

"""
    hyperbolic_form(kk)
    hyperbolic_form(kk, n)

The [`GWClass`](@ref) of the hyperbolic form ``\\mathbb{H} = ⟨1, -1⟩`` over the field
or finite étale algebra `kk`, or of the totally hyperbolic form
``(n/2)·\\mathbb{H}`` when an (even) rank `n` is specified. Odd `n` is an error.

# Examples
```julia-repl
julia> hyperbolic_form(GF(7))
[1   0]
[0   6]

julia> hyperbolic_form(Float64, 4)
4×4 Matrix{Float64}:
 1.0   0.0  0.0   0.0
 0.0  -1.0  0.0   0.0
 0.0   0.0  1.0   0.0
 0.0   0.0  0.0  -1.0
```

See also [`diagonal_form`](@ref), [`hyperbolic_unstable_form`](@ref),
[`sum_decomposition`](@ref).
"""
hyperbolic_form(kk) = hyperbolic_form(kk, 2)
function hyperbolic_form(kk, n::Integer)
    isodd(n) && error("entered rank is odd")
    _diagonal_form(kk, [isodd(i) ? 1 : -1 for i in 1:n])
end

# --- M2: makePfisterForm ----------------------------------------------------------
# <<a_1,...,a_n>> = <1,-a_1> ⊗ ... ⊗ <1,-a_n>

function _pfister_form(kk, entries)
    out = diagonal_form(kk, 1)
    for b in entries
        c = _coerce_into_or_error(kk, b)
        out = gw_tensor_product(out, diagonal_form(kk, (1, -c)))
    end
    out
end

"""
    pfister_form(kk, (a₁, …, aₙ))
    pfister_form(kk, a₁, a₂, …)

The [`GWClass`](@ref) of the Pfister form ``⟨⟨a_1, …, a_n⟩⟩`` over the field
`kk`: the rank-``2^n`` tensor product
``⟨1, -a_1⟩ ⊗ ⋯ ⊗ ⟨1, -a_n⟩``. A single entry produces a one-fold Pfister
form.

# Examples
```julia-repl
julia> pfister_form(QQ, (2, 6))
[1    0    0    0]
[0   -6    0    0]
[0    0   -2    0]
[0    0    0   12]

julia> pfister_form(GF(13), -2//3)
[1   0]
[0   5]
```

See also [`diagonal_form`](@ref), [`gw_tensor_product`](@ref).
"""
pfister_form(kk, L::Tuple) = _pfister_form(kk, collect(L))
pfister_form(kk, entries...) = _pfister_form(kk, collect(entries))

# --- M2: makeDiagonalUnstableForm / makeHyperbolicUnstableForm ---------------------

"""
    diagonal_unstable_form(kk, (a₁, …, aₙ))
    diagonal_unstable_form(kk, a₁, a₂, …)

The [`GWuClass`](@ref) represented by the diagonal form ``⟨a_1, …, a_n⟩``
over the field or finite étale algebra `kk`, with scalar the determinant
``a_1 ⋯ a_n``. See [`diagonal_form`](@ref) for the stable counterpart.

# Examples
```julia-repl
julia> diagonal_unstable_form(QQ, (3, 5, 7))
([3 0 0; 0 5 0; 0 0 7], 105)
```
"""
diagonal_unstable_form(kk, L::Tuple) = GWuClass(diagonal_form(kk, L))
diagonal_unstable_form(kk, entries...) = GWuClass(_diagonal_form(kk, collect(entries)))

"""
    hyperbolic_unstable_form(kk)
    hyperbolic_unstable_form(kk, n)

The [`GWuClass`](@ref) represented by the hyperbolic form ``\\mathbb{H} = ⟨1, -1⟩``
(or the totally hyperbolic form of even rank `n`) over the field or finite
étale algebra `kk`, with scalar its determinant. See
[`hyperbolic_form`](@ref) for the stable counterpart.

# Examples
```julia-repl
julia> hyperbolic_unstable_form(GF(7))
([1 0; 0 6], 6)
```
"""
hyperbolic_unstable_form(kk) = GWuClass(hyperbolic_form(kk))
hyperbolic_unstable_form(kk, n::Integer) = GWuClass(hyperbolic_form(kk, n))
