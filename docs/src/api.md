```@meta
CurrentModule = MotivicHomotopy
```

# API reference

All exported functions, grouped by task. Only user-facing functions are
listed; backend helpers are internal.

```@contents
Pages = ["api.md"]
Depth = 2
```

## Grothendieck–Witt classes

Constructing stable and unstable classes and combining them with the
Grothendieck–Witt ring operations.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["gw_classes.jl", "gwu_classes.jl"]
```

## Working with forms

Building standard forms, diagonalizing them, extracting simplified
representatives, and computing the Witt (sum) decomposition and anisotropic
part.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["building_forms.jl", "matrix_methods.jl", "simplified_representatives.jl", "decomposition.jl"]
```

## Invariants and classification

Numerical and arithmetic invariants of a form (rank, signature, discriminant,
Hilbert symbols, Hasse–Witt invariants, anisotropic dimension, Witt index),
the isotropy predicates, and the isomorphism test that classifies forms up to
equivalence.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["invariants.jl", "hilbert_symbols.jl", "arithmetic.jl", "anisotropic_dimension.jl", "isotropy.jl", "isomorphism.jl"]
```

## Computing ``\mathbb{A}^1``-degrees

The local and global ``\mathbb{A}^1``-Brouwer degrees, both stable and unstable.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["degrees.jl", "unstable_degrees.jl"]
```

## Étale algebras and transfer

Multiplication matrices, trace and norm of an étale algebra, and the transfer
(corestriction) of a Grothendieck–Witt class.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["trace_norm.jl", "transfer.jl"]
```
