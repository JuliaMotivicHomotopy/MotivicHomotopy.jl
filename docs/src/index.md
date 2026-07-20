```@meta
CurrentModule = MotivicHomotopy
```

# MotivicHomotopy.jl

```@docs
MotivicHomotopy
```

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl")
```

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

The *unstable* degrees compute numerically over ``\mathbb{C}`` when given
`HomotopyContinuation` expressions (the stable degrees are exact everywhere —
over ``\mathbb{C}`` they are just the identity form of the algebra dimension):

```julia
import HomotopyContinuation; HomotopyContinuation.@var x
global_unstable_A1_degree((x-1)*(x-2)*(x-3), (x-1)*(x-4))
```

Over ``\mathbb{R}`` the degree functions ask you to compute over ``\mathbb{Q}`` and base-change the
result. Forms over ``\mathbb{R}`` / ``\mathbb{C}`` are represented by `Matrix{Float64}` /
`Matrix{ComplexF64}` Gram matrices; pass `Float64` / `ComplexF64` where a
base field is expected.

See the [API reference](@ref) for every exported function, grouped by topic.

## Conventions

- **`==` is literal Gram-matrix equality**; mathematical equality in
  ``\text{GW}(k)`` is [`is_isomorphic_form`](@ref).
- Accessors and operations whose natural names would collide with Oscar
  exports carry a `gw_` prefix ([`gw_matrix`](@ref), not `gram_matrix`;
  [`gw_witt_index`](@ref), not `witt_index`).
- Degree Gram matrices are expressed in a standard-monomial basis; an
  equivalent form in another basis is the same class, so compare with
  [`is_isomorphic_form`](@ref) rather than entrywise.

## Version history

MotivicHomotopy.jl began as a Julia port of the Macaulay2 package
A1BrouwerDegrees. Its initial release consolidates the functionality of that
package's versions 1.1 and 2.0 into a single starting point; subsequent
releases listed below build on it.

- **v0.1.0** — based on versions 1.1 and 2.0 of the Macaulay2 package
  A1BrouwerDegrees. Implements the computation of local and global 𝔸¹-Brouwer
  degrees and of Grothendieck–Witt classes and their invariants; the
  computation of unstable local and global 𝔸¹-Brouwer degrees and manipulation
  of the unstable Grothendieck–Witt group; and the extension of
  Grothendieck–Witt class manipulations over fields to finite étale algebras
  over fields, including transfers along finite étale extensions. Developed by
  S. Atherton, N. Borisov, T. Brazelton, S. Dutta, F. Espino, T. Hagedorn,
  Z. Han, J. Lopez Garcia, J. Louwsma, Y. Luo, G. Ong, R. Sagayaraj, and
  A. Tawfeek.
