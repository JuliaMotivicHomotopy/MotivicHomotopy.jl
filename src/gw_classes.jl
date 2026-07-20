# Port of M2 Code/GrothendieckWittClasses.m2.
#
# Representation: M2's GrothendieckWittClass is a HashTable {matrix, cache}.
# Here: a struct holding the Gram matrix plus a Dict cache (M2 caches only the
# diagonal class; the Dict reproduces that mutable cache on an immutable class).
#
# M2's RR_53 / CC_53 are represented by Matrix{Float64} / Matrix{ComplexF64};
# exact fields and étale algebras use Oscar matrices (MatElem).

const FloatGram = Union{Matrix{Float64}, Matrix{ComplexF64}}

"""
    GWClass(M)

The isomorphism class of the nondegenerate symmetric bilinear form with Gram
matrix `M`, as an element of the Grothendieck–Witt ring ``\\text{GW}(k)`` of a field
(or finite étale algebra over a field) of characteristic not 2.

Given a basis ``e_1, …, e_n`` of a ``k``-vector space ``V``, a symmetric
bilinear form ``β : V × V → k`` is encoded by its Gram matrix
``(β(e_i, e_j))_{i,j}``; a change of basis replaces the Gram matrix by a
congruent one, so a symmetric matrix determines the form up to congruence.
The constructor checks that `M` is symmetric, nondegenerate, and defined over
a supported coefficient ring, and errors otherwise.

Supported Gram matrices:
- an Oscar matrix (`MatElem`) over `QQ`, a finite field of odd characteristic,
  or a finite étale algebra (a zero-dimensional `MPolyQuoRing` with
  nondegenerate trace form);
- a plain Julia matrix of real or complex numbers, stored as `Matrix{Float64}`
  or `Matrix{ComplexF64}` — these play the role of forms over ``\\mathbb{R}`` and ``\\mathbb{C}``.

Equality `==` compares Gram matrices literally (same base ring, same entries).
Mathematical equality in ``\\text{GW}(k)`` is tested with
[`is_isomorphic_form`](@ref).

The Gram matrix, coefficient algebra, and base field are recovered with
[`gw_matrix`](@ref), [`gw_algebra`](@ref), and [`gw_base_field`](@ref); a
diagonal representative with [`diagonal_class`](@ref). Further invariants:
[`form_rank`](@ref), [`form_signature`](@ref),
[`integral_discriminant`](@ref), [`hasse_witt_invariant`](@ref),
[`anisotropic_dimension`](@ref), [`anisotropic_part`](@ref),
[`sum_decomposition`](@ref), [`is_isotropic_form`](@ref),
[`is_anisotropic_form`](@ref).

# Examples
```julia-repl
julia> beta = GWClass(QQ[2 1; 1 3])
[2   1]
[1   3]

julia> gw_matrix(beta)
[2   1]
[1   3]

julia> gw_base_field(beta)
Rational field

julia> GWClass([0.0 1.0; 1.0 0.0])    # a form over the real numbers
2×2 Matrix{Float64}:
 0.0  1.0
 1.0  0.0
```
"""
struct GWClass{T}
    matrix::T
    cache::Dict{Symbol, Any}
end

# --- well-definedness (M2: isWellDefinedGW) ------------------------------------

_ring_char(R) = characteristic(R)
_ring_char(R::MPolyQuoRing) = characteristic(coefficient_ring(base_ring(R)))

function _is_well_defined_gw(M::MatElem)
    _is_symmetric_square(M) || return false
    iszero(det(M)) && return false
    R = base_ring(M)
    (R isa Field || is_finite_etale_algebra(R)) || return false
    _ring_char(R) != 2 || return false
    true
end

function _is_well_defined_gw(M::FloatGram)
    _is_symmetric_square(M) || return false
    !iszero(LinearAlgebra.det(M))
end

const _GW_ERR = "makeGWClass called on a matrix that does not represent a " *
    "nondegenerate symmetric bilinear form over a field of characteristic not 2"

# --- constructors (M2: makeGWClass) --------------------------------------------

function GWClass(M::MatElem)
    _is_well_defined_gw(M) || error(_GW_ERR)
    GWClass{typeof(M)}(M, Dict{Symbol, Any}())
end

function GWClass(M::AbstractMatrix{<:Real})
    Mf = Matrix{Float64}(M)
    _is_well_defined_gw(Mf) || error(_GW_ERR)
    GWClass{Matrix{Float64}}(Mf, Dict{Symbol, Any}())
end

function GWClass(M::AbstractMatrix{<:Complex})
    Mf = Matrix{ComplexF64}(M)
    _is_well_defined_gw(Mf) || error(_GW_ERR)
    GWClass{Matrix{ComplexF64}}(Mf, Dict{Symbol, Any}())
end

# --- accessors (M2: getMatrix / getAlgebra / getBaseField) ----------------------

"""
    gw_matrix(beta)

The Gram matrix of a [`GWClass`](@ref) or [`GWuClass`](@ref): a symmetric
matrix over the coefficient algebra of the class (`MatElem` for exact fields
and étale algebras, `Matrix{Float64}` / `Matrix{ComplexF64}` over ``\\mathbb{R}`` / ``\\mathbb{C}``).

# Examples
```julia-repl
julia> beta = GWClass(QQ[2 1; 1 3]);

julia> gw_matrix(beta)
[2   1]
[1   3]
```

See also [`gw_scalar`](@ref), [`gw_algebra`](@ref), [`gw_base_field`](@ref).
"""
gw_matrix(beta::GWClass) = beta.matrix

_algebra_of(M::MatElem) = base_ring(M)
_algebra_of(::Matrix{Float64}) = Float64
_algebra_of(::Matrix{ComplexF64}) = ComplexF64

function _base_field_of(M::MatElem)
    R = base_ring(M)
    (R isa QQField || R isa FinField) && return R
    R isa Field && return R
    if R isa MPolyQuoRing
        # M2: requires the zero ideal to be prime, then applies toField.
        # Julia has no toField; we return the (domain) quotient ring itself.
        is_prime(ideal(R, [zero(R)])) ||
            error("the Grothendieck-Witt class is not defined over a field")
        return R
    end
    error("the Grothendieck-Witt class is not defined over a field")
end
_base_field_of(::Matrix{Float64}) = Float64
_base_field_of(::Matrix{ComplexF64}) = ComplexF64

"""
    gw_algebra(beta)

The coefficient algebra over which a [`GWClass`](@ref) or [`GWuClass`](@ref)
is defined — a field or a finite étale algebra over a field. Over ``\\mathbb{R}`` / ``\\mathbb{C}``
(float-backed classes) this returns the type `Float64` / `ComplexF64`.

# Examples
```julia-repl
julia> gw_algebra(GWClass(QQ[2 1; 1 3]))
Rational field
```

See also [`gw_base_field`](@ref), [`gw_matrix`](@ref).
"""
gw_algebra(beta::GWClass) = _algebra_of(beta.matrix)

"""
    gw_base_field(beta)

The base field of a [`GWClass`](@ref) or [`GWuClass`](@ref), when the class is
defined over a field. For a class over an étale algebra this checks that the
zero ideal is prime (i.e. the algebra is a field) and errors otherwise; for a
class that already lives over `QQ` or a finite field it returns that field.
Over ``\\mathbb{R}`` / ``\\mathbb{C}`` (float-backed classes) it returns the type `Float64` /
`ComplexF64`.

# Examples
```julia-repl
julia> gw_base_field(GWClass(QQ[2 1; 1 3]))
Rational field
```

See also [`gw_algebra`](@ref).
"""
gw_base_field(beta::GWClass) = _base_field_of(beta.matrix)

# --- equality (M2 `===`: literal Gram matrix, cache ignored) --------------------

_gram_equal(A::MatElem, B::MatElem) =
    base_ring(A) === base_ring(B) && size(A) == size(B) && A == B
_gram_equal(A::Matrix{T}, B::Matrix{T}) where {T <: Union{Float64, ComplexF64}} =
    size(A) == size(B) && A == B
_gram_equal(::Any, ::Any) = false

Base.:(==)(a::GWClass, b::GWClass) = _gram_equal(a.matrix, b.matrix)
Base.hash(a::GWClass, h::UInt) = hash(a.matrix, hash(:GWClass, h))

# M2 prints the underlying matrix (net GrothendieckWittClass := net getMatrix)
Base.show(io::IO, beta::GWClass) = show(io, beta.matrix)
Base.show(io::IO, mime::MIME"text/plain", beta::GWClass) =
    show(io, mime, beta.matrix)

# --- operations (M2: addGW / multiplyGW) ----------------------------------------

const _FIELD_MISMATCH_ERR = "these classes have different underlying fields"

# M2 allows two Galois fields of the same order; map the second into the first.
function _align_second(Mb::MatElem, Mg::MatElem)
    Kb, Kg = base_ring(Mb), base_ring(Mg)
    if Kb isa FinField && Kg isa FinField
        order(Kb) == order(Kg) || error(_FIELD_MISMATCH_ERR)
        Kb === Kg && return Mg
        # prime fields transfer through ZZ (M2: substitute(-, Kb)); a same-order
        # transfer between distinct non-prime fields would need an embedding.
        return map_entries(x -> Kb(lift(ZZ, x)), Mg)
    end
    Kb === Kg || error(_FIELD_MISMATCH_ERR)
    Mg
end

"""
    gw_direct_sum(beta, gamma)

The direct (block) sum of two classes over the same base field — addition in
the Grothendieck–Witt ring. For two [`GWClass`](@ref) inputs the result has
the block-diagonal Gram matrix; for two [`GWuClass`](@ref) inputs the Gram
matrices are block-summed and the scalars multiplied. `beta + gamma` is an
alias.

# Examples
```julia-repl
julia> beta1 = GWClass(QQ[1 2; 2 3]);

julia> beta2 = GWClass(QQ[3 4; 4 5]);

julia> gw_direct_sum(beta1, beta2)
[1   2   0   0]
[2   3   0   0]
[0   0   3   4]
[0   0   4   5]

julia> gw_direct_sum(GWuClass(QQ[2 1; 1 2]), GWuClass(QQ[1 2; 2 6]))
([2 1 0 0; 1 2 0 0; 0 0 1 2; 0 0 2 6], 6)
```

See also [`gw_tensor_product`](@ref), [`divisorial_sum`](@ref).
"""
function gw_direct_sum(beta::GWClass{<:MatElem}, gamma::GWClass{<:MatElem})
    Mg = _align_second(beta.matrix, gamma.matrix)
    GWClass(block_diagonal_matrix([beta.matrix, Mg]))
end

"""
    gw_tensor_product(beta, gamma)

The tensor product of two [`GWClass`](@ref)es over the same base field —
multiplication in the Grothendieck–Witt ring. The resulting Gram matrix is
the Kronecker product of the two Gram matrices. `beta * gamma` is an alias.

# Examples
```julia-repl
julia> beta1 = GWClass(QQ[1 2; 2 3]);

julia> beta2 = GWClass(QQ[3 4; 4 5]);

julia> gw_tensor_product(beta1, beta2)
[3    4    6    8]
[4    5    8   10]
[6    8    9   12]
[8   10   12   15]
```

See also [`gw_direct_sum`](@ref).
"""
function gw_tensor_product(beta::GWClass{<:MatElem}, gamma::GWClass{<:MatElem})
    Mg = _align_second(beta.matrix, gamma.matrix)
    GWClass(kronecker_product(beta.matrix, Mg))
end

function _float_block_diag(A::Matrix{T}, B::Matrix{T}) where {T}
    n1, n2 = size(A, 1), size(B, 1)
    C = zeros(T, n1 + n2, n1 + n2)
    C[1:n1, 1:n1] = A
    C[(n1 + 1):end, (n1 + 1):end] = B
    C
end

gw_direct_sum(beta::GWClass{Matrix{T}}, gamma::GWClass{Matrix{T}}) where {T} =
    GWClass(_float_block_diag(beta.matrix, gamma.matrix))
gw_tensor_product(beta::GWClass{Matrix{T}}, gamma::GWClass{Matrix{T}}) where {T} =
    GWClass(LinearAlgebra.kron(beta.matrix, gamma.matrix))

gw_direct_sum(::GWClass, ::GWClass) = error(_FIELD_MISMATCH_ERR)
gw_tensor_product(::GWClass, ::GWClass) = error(_FIELD_MISMATCH_ERR)

# Optional operator aliases; the named functions are canonical.
Base.:+(beta::GWClass, gamma::GWClass) = gw_direct_sum(beta, gamma)
Base.:*(beta::GWClass, gamma::GWClass) = gw_tensor_product(beta, gamma)
