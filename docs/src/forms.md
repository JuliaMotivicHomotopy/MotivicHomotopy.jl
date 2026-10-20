```@meta
CurrentModule = MotivicHomotopy
```

# Working with forms

Everything for building and computing with Grothendieck–Witt classes:
constructors, the ring operations, simplified representatives, the Witt
decomposition and anisotropic part, invariants, and the isomorphism test.

Supported base fields: ``\mathbb{Q}``, finite fields of odd characteristic,
and ``\mathbb{R}`` / ``\mathbb{C}`` via `Matrix{Float64}` /
`Matrix{ComplexF64}` Gram matrices (pass `Float64` / `ComplexF64` where a
base field is expected).

## Classes and ring operations

Constructing stable and unstable classes, their accessors, and combining
them with the Grothendieck–Witt ring operations (orthogonal sum and tensor
product) and the divisorial sum of unstable local degrees.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["gw_classes.jl", "gwu_classes.jl"]
```

## Building and decomposing forms

Building standard forms, diagonalizing them by congruence, extracting
simplified representatives, and computing the Witt (sum) decomposition and
anisotropic part.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["building_forms.jl", "matrix_methods.jl", "simplified_representatives.jl", "decomposition.jl"]
```

## Invariants and classification

Rank, signature, discriminant, Hasse–Witt invariants, anisotropic dimension
and Witt index, the isotropy predicates, and the isomorphism test that
classifies forms up to equivalence. `==` is literal Gram-matrix equality;
mathematical equality in ``\text{GW}(k)`` is tested by
[`is_isomorphic_form`](@ref).

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["invariants.jl", "anisotropic_dimension.jl", "isotropy.jl", "isomorphism.jl"]
```

## Étale algebras and transfer

Multiplication matrices, trace and norm of an étale algebra, and the
transfer (corestriction) of a Grothendieck–Witt class.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["trace_norm.jl", "transfer.jl"]
```
