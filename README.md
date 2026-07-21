<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo-dark.png">
    <img src="docs/src/assets/logo.png" alt="MotivicHomotopy.jl" width="480">
  </picture>
</p>

<p align="center">
  <a href="https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"></a>
  <a href="https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"></a>
  <a href="https://github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"></a>
  <a href="https://codecov.io/gh/JuliaMotivicHomotopy/MotivicHomotopy.jl"><img src="https://codecov.io/gh/JuliaMotivicHomotopy/MotivicHomotopy.jl/branch/main/graph/badge.svg" alt="Coverage"></a>
</p>

A Julia package for computing local and global $`\mathbb{A}^1`$-Brouwer degrees
and studying the symmetric bilinear forms they produce, over $`\mathbb{Q}`$,
finite fields of odd characteristic, $`\mathbb{R}`$, and $`\mathbb{C}`$.

The $`\mathbb{A}^1`$-Brouwer degree of an endomorphism of affine space with
isolated zeros is a class in the Grothendieck–Witt ring $`\text{GW}(k)`$ of
symmetric bilinear forms: its rank recovers the degree of the associated
complex map and its signature the degree of the associated real map. For
pointed rational functions $`\mathbb{P}^1 \to \mathbb{P}^1`$ the degree lives in
the unstable group $`\text{GW}^u(k)`$, which additionally records a
$`k^\times`$-scalar.

## Features

- **Grothendieck–Witt classes** — `GWClass` and unstable `GWuClass` with
  accessors, direct sum, tensor product, and the divisorial sum of unstable
  local degrees.
- **Form constructors** — diagonal, hyperbolic, and Pfister forms, with
  unstable variants.
- **Simplification** — diagonalization via congruence, square-class-reduced
  diagonal representatives, Witt (sum) decomposition, and the anisotropic part.
- **Invariants and classification** — rank, signature, integral discriminant,
  Hilbert symbols, Hasse–Witt invariants, anisotropic dimension, Witt index,
  isotropy predicates, and the isomorphism test `is_isomorphic_form`.
- **$`\mathbb{A}^1`$-Brouwer degrees** — local and global, stable and unstable,
  computed via the multivariate Bézoutian bilinear form.
- **Finite étale algebras** — multiplication matrices, trace and norm, and the
  transfer map $`\text{GW}(L) \to \text{GW}(k)`$ along a finite étale extension.

## Installation

The package requires Julia ≥ 1.11 and is not yet registered; install it
directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl")
```

Its main computational dependencies, [Oscar](https://www.oscar-system.org/) and
[HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/), are
installed automatically.

## Quick start

```julia
using Oscar, MotivicHomotopy

# The A¹-degree of z ↦ z² : rank 2 (complex degree), signature 0 (real degree)
S, (x,) = polynomial_ring(QQ, ["x"])
beta = global_A1_degree([x^2 + 1])
form_rank(beta), form_signature(beta)        # (2, 0)

# Local degrees sum to the global degree
S, (x, y) = polynomial_ring(QQ, ["x", "y"])
f = [x^3 - x^2 - y, y]
d1 = local_A1_degree(f, ideal(S, [x - 1, y]))
d2 = local_A1_degree(f, ideal(S, [x, y]))
is_isomorphic_form(global_A1_degree(f), gw_direct_sum(d1, d2))   # true

# Witt decomposition of a form
sum_decomposition_string(GWClass(QQ[1 2 3; 2 4 5; 3 5 6]))       # "H + <1>"
```

The unstable degrees compute numerically over $`\mathbb{C}`$ when given
`HomotopyContinuation` expressions:

```julia
import HomotopyContinuation; HomotopyContinuation.@var x
global_unstable_A1_degree((x-1)*(x-2)*(x-3), (x-1)*(x-4))
```

## Base fields

Which arithmetic runs is determined by the base field:

| Base field | Computation |
|---|---|
| $`\mathbb{Q}`$ and $`\mathbb{F}_q`$ (odd $`q`$) | exact, via Oscar |
| $`\mathbb{R}`$ | compute over $`\mathbb{Q}`$ and base-change the resulting Gram matrix (`Matrix{Float64}` represents forms over $`\mathbb{R}`$) |
| $`\mathbb{C}`$ | stable degrees are rank-only and come from the exact computation; unstable degrees are computed numerically via HomotopyContinuation.jl |

See the [documentation](https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl/dev/)
for the full API reference, conventions (e.g. `==` is literal Gram-matrix
equality; mathematical equality is `is_isomorphic_form`), and version history.

## Relation to the Macaulay2 package

MotivicHomotopy.jl began as a Julia port of the Macaulay2 package
[A1BrouwerDegrees](https://msp.org/jsag/2024/14-1/p15.xhtml) (v1.1 published in
the Journal of Software for Algebra and Geometry 14, 2024); its initial release
consolidates the functionality of that package's versions 1.1 and 2.0 into a
single starting point. The Julia package has its own version lineage — see
[CONTRIBUTING.md](CONTRIBUTING.md) for the versioning policy.

## Background

- F. Morel, *$`\mathbb{A}^1`$-algebraic topology over a field*, Springer Lecture Notes in Mathematics, 2012.
- C. Cazanave, *Algebraic homotopy classes of rational functions*, Annales Scientifiques de l'École Normale Supérieure, 2012.
- J. L. Kass and K. Wickelgren, *The class of Eisenbud–Khimshiashvili–Levine is the local $`\mathbb{A}^1`$-Brouwer degree*, Duke Mathematical Journal, 2019.
- T. Brazelton, S. McKean, and S. Pauli, *Bézoutians and the $`\mathbb{A}^1`$-degree*, Algebra & Number Theory, 2023.

## Acknowledgments

This package began as a Julia port of the `A1BrouwerDegrees` Macaulay2 package
(v1 + v2). The port was carried out with the assistance of Claude Code.
