# Port of M2 Code/UnstableGrothendieckWittClasses.m2 (plus isCompatibleElement
# from Code/BuildingForms.m2, folded into scalar coercion).
#
# Representation: M2's UnstableGrothendieckWittClass is a HashTable
# {matrix, cache, scalar}; here a struct with Gram matrix, scalar (an element of
# the same ring), and the Dict cache.

"""
    GWuClass(M)
    GWuClass(M, a)
    GWuClass(beta)
    GWuClass(beta, a)

An element of the unstable Grothendieck–Witt group
``\\text{GW}^u(k) = \\text{GW}(k) ×_{k^×/(k^×)^2} k^×`` of a field (or finite étale algebra)
of characteristic not 2: the data of a [`GWClass`](@ref) together with a
nonzero scalar whose square class agrees with the determinant of the Gram
matrix.

The constructor takes a symmetric matrix `M` (or an existing `GWClass`
`beta`) and optionally a scalar `a`; when the scalar is omitted, the
determinant of the Gram matrix is used. Over ``\\mathbb{Q}``, ``\\mathbb{R}``, and finite fields of odd
characteristic the constructor verifies that `a` agrees with the determinant
up to squares and errors otherwise. Over an arbitrary finite étale algebra
this verification is not possible; any nonzero scalar is accepted and a
warning is printed when it differs from the determinant — the user must check
the square-class condition by hand.

Accessors: [`gw_matrix`](@ref), [`gw_scalar`](@ref), [`gw_algebra`](@ref),
[`gw_base_field`](@ref), and [`stable_part`](@ref) for the underlying
``\\text{GW}(k)``-class. Operations: [`gw_direct_sum`](@ref) (alias `+`) and
[`divisorial_sum`](@ref). Equality `==` compares the Gram matrix and the
scalar literally; use [`is_isomorphic_form`](@ref) for equality in
``\\text{GW}^u(k)``.

# Examples
```julia-repl
julia> M = QQ[0 1; 1 0];

julia> alpha = GWuClass(M, -4)
([0 1; 1 0], -4)

julia> gw_scalar(alpha)
-4

julia> GWuClass(M)               # scalar defaults to det(M)
([0 1; 1 0], -1)

julia> GWuClass(GWClass(M), -9)
([0 1; 1 0], -9)
```
"""
struct GWuClass{T, S}
    matrix::T
    scalar::S
    cache::Dict{Symbol, Any}
end

# --- scalar compatibility + coercion (M2: isCompatibleElement) -------------------
# Returns the scalar as an element of the matrix's ring, or `nothing` if the
# scalar is not compatible (M2 then fails the well-definedness check).
# The ring-level coercion `_coerce_into` lives with the form builders
# (building_forms.jl), matching M2's placement of isCompatibleElement.

_coerce_scalar(M::MatElem, b) = _coerce_into(base_ring(M), b)
_coerce_scalar(::Matrix{Float64}, b) = _coerce_into(Float64, b)
_coerce_scalar(::Matrix{ComplexF64}, b) = _coerce_into(ComplexF64, b)

# --- well-definedness (M2: isWellDefinedGWu) -------------------------------------

# unit test: fields need nonzero; étale algebras need 1 ∈ I + (c)
_is_unit_scalar(c::FieldElem) = !iszero(c)
function _is_unit_scalar(c::MPolyQuoRingElem)
    A = parent(c)
    R = base_ring(A)
    is_one(ideal(R, vcat(gens(modulus(A)), [lift(c)])))
end
_is_unit_scalar(c::Union{Float64, ComplexF64}) = !iszero(c)

function _is_well_defined_gwu(M::MatElem, c)
    _is_unit_scalar(c) || return false
    R = base_ring(M)
    d = det(M)
    if R isa QQField
        # M2 compares squarefree parts; equivalent: d and c are in the same square
        # class over QQ iff d*c is a square.
        is_square(d * c) || return false
    elseif R isa FinField
        is_square(d) == is_square(c) || return false
    else
        # étale algebra: M2 cannot verify the square class; it warns when c != det.
        d == c || @warn "Warning, unable to verify whether the determinant of M " *
            "and b agree up to squares."
    end
    _is_well_defined_gw(M)
end

function _is_well_defined_gwu(M::Matrix{Float64}, c::Float64)
    iszero(c) && return false
    sign(c) == sign(LinearAlgebra.det(M)) || return false   # M2: signs must agree over RR
    _is_well_defined_gw(M)
end

function _is_well_defined_gwu(M::Matrix{ComplexF64}, c::ComplexF64)
    iszero(c) && return false                               # M2: scalar must be a unit
    _is_well_defined_gw(M)
end

const _GWU_ERR = "makeGWuClass called on a pair that does not produce a " *
    "well-defined element of the unstable Grothendieck-Witt group."

function _make_gwu(M, b)
    c = _coerce_scalar(M, b)
    c === nothing && error(_GWU_ERR)
    _is_well_defined_gwu(M, c) || error(_GWU_ERR)
    GWuClass{typeof(M), typeof(c)}(M, c, Dict{Symbol, Any}())
end

# --- constructors (M2: makeGWuClass, six variants) -------------------------------

_gram_det(M::MatElem) = det(M)
_gram_det(M::FloatGram) = LinearAlgebra.det(M)

GWuClass(M::MatElem) = _make_gwu(M, det(M))
GWuClass(M::MatElem, b) = _make_gwu(M, b)
GWuClass(M::AbstractMatrix{<:Real}) =
    (Mf = Matrix{Float64}(M); _make_gwu(Mf, LinearAlgebra.det(Mf)))
GWuClass(M::AbstractMatrix{<:Real}, b) = _make_gwu(Matrix{Float64}(M), b)
GWuClass(M::AbstractMatrix{<:Complex}) =
    (Mf = Matrix{ComplexF64}(M); _make_gwu(Mf, LinearAlgebra.det(Mf)))
GWuClass(M::AbstractMatrix{<:Complex}, b) = _make_gwu(Matrix{ComplexF64}(M), b)
GWuClass(beta::GWClass) = _make_gwu(beta.matrix, _gram_det(beta.matrix))
GWuClass(beta::GWClass, b) = _make_gwu(beta.matrix, b)

# --- accessors (M2: getMatrix / getScalar / getAlgebra / getBaseField / getGWClass)

gw_matrix(beta::GWuClass) = beta.matrix

"""
    gw_scalar(beta::GWuClass)

The ``k^×``-factor of an unstable Grothendieck–Witt class: the nonzero
scalar that, together with the Gram matrix, determines the class. It is an
element of the coefficient algebra of `beta`.

# Examples
```julia-repl
julia> gw_scalar(GWuClass(QQ[0 1; 1 0], -4))
-4
```

See also [`GWuClass`](@ref), [`gw_matrix`](@ref), [`stable_part`](@ref).
"""
gw_scalar(beta::GWuClass) = beta.scalar

gw_algebra(beta::GWuClass) = _algebra_of(beta.matrix)
gw_base_field(beta::GWuClass) = _base_field_of(beta.matrix)

"""
    stable_part(beta::GWuClass)

The image of an unstable Grothendieck–Witt class under the projection
``\\text{GW}^u(k) → \\text{GW}(k)``: the [`GWClass`](@ref) with the same Gram matrix,
forgetting the scalar.

# Examples
```julia-repl
julia> stable_part(GWuClass(QQ[0 1; 1 0], -4))
[0   1]
[1   0]
```

See also [`GWuClass`](@ref), [`gw_scalar`](@ref).
"""
stable_part(beta::GWuClass) = GWClass(beta.matrix)

# --- equality (M2 `===`: matrix and scalar, cache ignored) ------------------------

_scalar_equal(x::RingElem, y::RingElem) = parent(x) === parent(y) && x == y
_scalar_equal(x::T, y::T) where {T <: Union{Float64, ComplexF64}} = x == y
_scalar_equal(::Any, ::Any) = false

Base.:(==)(a::GWuClass, b::GWuClass) =
    _gram_equal(a.matrix, b.matrix) && _scalar_equal(a.scalar, b.scalar)
Base.hash(a::GWuClass, h::UInt) = hash(a.scalar, hash(a.matrix, hash(:GWuClass, h)))

# M2 prints the (matrix, scalar) pair
Base.show(io::IO, beta::GWuClass) = print(io, "(", beta.matrix, ", ", beta.scalar, ")")

# --- addGWu -----------------------------------------------------------------------
# NOTE: M2's addGWu compares getBaseField (not getAlgebra), so it errors over
# non-domain étale algebras; ported faithfully.

function gw_direct_sum(beta::GWuClass{<:MatElem}, gamma::GWuClass{<:MatElem})
    gw_base_field(beta); gw_base_field(gamma)   # M2 evaluates these (may error)
    Mg = _align_second(beta.matrix, gamma.matrix)
    sg = base_ring(Mg) === parent(gamma.scalar) ? gamma.scalar :
        base_ring(Mg)(lift(ZZ, gamma.scalar))
    GWuClass(block_diagonal_matrix([beta.matrix, Mg]), beta.scalar * sg)
end

gw_direct_sum(beta::GWuClass{Matrix{T}}, gamma::GWuClass{Matrix{T}}) where {T} =
    GWuClass(_float_block_diag(beta.matrix, gamma.matrix), beta.scalar * gamma.scalar)

gw_direct_sum(::GWuClass, ::GWuClass) = error(_FIELD_MISMATCH_ERR)

Base.:+(beta::GWuClass, gamma::GWuClass) = gw_direct_sum(beta, gamma)

# --- addGWuDivisorial -------------------------------------------------------------

_gram_rank(M::MatElem) = rank(M)
_gram_rank(M::FloatGram) = LinearAlgebra.rank(M)

"""
    divisorial_sum(class_list, root_list)

The divisorial sum of a list of unstable local degrees with respect to the
divisor of roots at which they were computed.

For a pointed rational function ``f/g : \\mathbb{P}^1_k → \\mathbb{P}^1_k`` with zeros
``r_1, …, r_n`` and unstable local ``\\mathbb{A}^1``-degrees ``β_1, …, β_n`` at those
zeros, the global unstable degree is *not* the [`gw_direct_sum`](@ref) of the
local degrees: the local-to-global formula of Igieobo et al. [I+24] weights the sum
by the configuration of the zeros. Concretely, the Gram matrices are
block-summed while the scalar picks up the factor
``∏_{i<j} (r_i - r_j)^{2 m_i m_j}`` (with ``m_i`` the rank of ``β_i``) on top
of the product of the scalars.

`class_list` is a vector of [`GWuClass`](@ref)es over a common base field and
`root_list` the corresponding vector of base-field elements.

# Examples

The local degrees of ``(x^2 + x - 2)/(3x + 5)`` at its zeros −2 and 1 are
``(⟨1/3⟩, 1/3)`` and ``(⟨8/3⟩, 8/3)``; their divisorial sum recovers the
global unstable degree:

```julia-repl
julia> alpha = GWuClass(matrix(QQ, 1, 1, [1//3]));

julia> beta = GWuClass(matrix(QQ, 1, 1, [8//3]));

julia> divisorial_sum([alpha, beta], [-2, 1])
([1//3 0; 0 8//3], 8)
```

# References
- [I+24] J. Igieobo et al., *Motivic configurations on the line*, Advances in Mathematics 482 (2025), 110637.

See also [`global_unstable_A1_degree`](@ref),
[`local_unstable_A1_degree`](@ref).
"""
function divisorial_sum(class_list::AbstractVector{<:GWuClass},
                        root_list::AbstractVector)
    n = length(class_list)
    n == length(root_list) || error("need same number of classes and roots")
    n == 0 && error("the empty sum is the additive identity of the unstable " *
        "Grothendieck-Witt group over the field of interest; please construct " *
        "this as GWuClass(matrix(k, 0, 0, []), 1)")

    # base fields must agree (Galois fields compare by order, as in M2)
    K1 = gw_base_field(class_list[1])
    for c in class_list
        K = gw_base_field(c)
        if K1 isa FinField
            (K isa FinField && order(K) == order(K1)) ||
                error("the list of GWu classes should have the same base field")
        else
            K === K1 || error("the list of GWu classes should have the same base field")
        end
    end

    # roots must lie in the base field (M2: isCompatibleElement per class)
    roots = Vector{Any}(undef, n)
    for i in 1:n
        r = _coerce_scalar(gw_matrix(class_list[i]), root_list[i])
        r === nothing && error("the roots must be in the base field of the classes")
        roots[i] = r
    end

    mults = [_gram_rank(gw_matrix(c)) for c in class_list]
    new_scalar = prod(gw_scalar(c) for c in class_list)
    for i in 1:n, j in (i + 1):n
        new_scalar *= (roots[i] - roots[j])^(2 * mults[i] * mults[j])
    end

    mats = [gw_matrix(c) for c in class_list]
    new_form = mats[1] isa MatElem ? block_diagonal_matrix(mats) :
        reduce(_float_block_diag, mats)
    GWuClass(new_form, new_scalar)
end
