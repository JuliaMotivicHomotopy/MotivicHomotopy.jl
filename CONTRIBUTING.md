# Contributing to MotivicHomotopy.jl

## Versioning

MotivicHomotopy.jl follows Semantic Versioning as interpreted by Julia's package
manager, and versions the package's public Julia API — not the maturity of the
underlying mathematics, which it inherits from the Macaulay2 package it is based
on.

The package is in the 0.x series deliberately: the API is still stabilizing.
Under Julia/Pkg's reading of SemVer, in 0.x the first nonzero component is the
breaking one. While the package is 0.MINOR.PATCH:

- 0.MINOR.PATCH → 0.MINOR.(PATCH+1): non-breaking — bug fixes and internal
  changes that do not alter the public API.
- 0.MINOR.PATCH → 0.(MINOR+1).0: may change the public API — new features
  and/or breaking changes.

Increments must be standard (the General registry requires it): bump exactly one
component by one, and never skip numbers. Any release that breaks the public API
must carry release notes that say so (mention "breaking").

1.0.0 is reserved for when the Julia API has stabilized — the core port is
complete, the type and naming conventions are settled, and the interface has
seen enough use to be worth committing to. Releasing 1.0.0 promises that
subsequent breaking changes will require 2.0.0.

### Relationship to A1BrouwerDegrees

MotivicHomotopy.jl has its own version lineage, independent of the Macaulay2
package A1BrouwerDegrees from which it originates. Its v0.1.0 consolidates the
functionality of that package's versions 1.1 and 2.0 as a single starting point;
this states provenance, not version inheritance. The M2 package's numbers are
not continued here.

### Keeping versions in sync

The version field in Project.toml is the single source of truth. Update the
"Version history" section of the documentation in the same change that bumps it,
so the two never drift.

## Naming conventions

To avoid clobbering or silently extending generic functions owned by Oscar or
Base, the package uses distinct, descriptive names for its own operations rather
than overloading existing generics. The rule is uniform: the package extends no
Oscar generic, and the only Base methods it defines are the standard interface
methods `==`, `hash`, and `show` plus the documented operator aliases `+` and
`*`, all on the package's own types (`GWClass` / `GWuClass`) — never on foreign
types.

Names used in place of generics that would otherwise collide include
`algebra_trace`, `algebra_norm`, `form_rank`, `form_signature`, and `gw_matrix`
(the Gram matrix of a Grothendieck–Witt class — named `gw_matrix` rather than
extending Oscar's `gram_matrix`). When adding functionality, prefer a distinct
name over a method on an Oscar/Base generic, even where extending the generic
would be type-safe.

## Documentation conventions

Mathematical notation follows the convention of the source package. In
docstrings, all mathematics is written as LaTeX inside double backticks;
Unicode glyphs (𝔸¹, ℚ, ...) are used only in code comments, not in docstrings.
In documentation pages, inline math uses double backticks and display math uses
fenced math blocks; the dollar-sign delimiters are not used.
