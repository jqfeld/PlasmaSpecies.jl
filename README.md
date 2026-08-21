# PlasmaSpecies

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jqfeld.github.io/PlasmaSpecies.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/PlasmaSpecies.jl/dev/)
[![Build Status](https://github.com/jqfeld/PlasmaSpecies.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jqfeld/PlasmaSpecies.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jqfeld/PlasmaSpecies.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jqfeld/PlasmaSpecies.jl)

## Overview

`PlasmaSpecies.jl` provides a small collection of types and helper functions for
working with the charged and neutral species that appear in plasma chemistry
models.  It focuses on

- convenient constructors for species such as electrons, ions, and neutrals using
  the [LoKI-B](https://github.com/IST-Lisbon/LoKI) string notation,
- utilities for navigating the hierarchy of vibrational, rotational, and
  electronic states via `SpeciesTree`, and
- parsing plasma reactions (including reversible ones) that can be converted into
  [Catalyst.jl](https://catalyst.sciml.ai/stable/) `Reaction`s.

Extensive API documentation and additional examples are available in the
[package documentation](https://jqfeld.github.io/PlasmaSpecies.jl/stable/).

## Installation

Install the package directly from the Git repository using the Julia package
manager:

```julia
julia> import Pkg; Pkg.add(url="https://github.com/jqfeld/PlasmaSpecies.jl.git")
```

The optional Catalyst.jl integration is automatically activated once
`using Catalyst` has been executed in the same Julia session.

## Quick start

```julia
julia> using PlasmaSpecies

julia> electron = p"e"
e

julia> nitrogen_ion = Species("N2[+,B]")
N2[+,B]

julia> gas(nitrogen_ion), charge(nitrogen_ion)
(N₂, +)

julia> reaction = p"e + N2 --> 2e + N2[+]"
e+N2-->2e+N2[+]

julia> reaction.subs, reaction.prods
([e, N2], [e, N2[+]])

julia> tree = SpeciesTree([electron, nitrogen_ion])
nothing
├─ e
└─ N2[+]
   └─ N2[+,B]
```

### Gases are compositions

A gas is resolved from its chemical formula, so anything built from known
nuclides works without being registered first — and knows its own mass:

```julia
julia> mass(Species("CO2")) / 1.66053906660e-27   # daltons
44.009

julia> gas(Species("H2O")), gas(Species("SF6"))
(H₂O, SF₆)
```

Write a specific isotopologue with ExoMol's `iso_slug` spelling, which lines
these strings up with `ExoMol.jl`'s `load_isotopologue`:

```julia
julia> mass(Species("28Si-16O")) / 1.66053906660e-27
43.97184115422001

julia> gas(Species("16O-1H"))
¹⁶O¹H
```

A bare element symbol means natural isotopic abundance (`N2` is 2 × 14.007 u),
which is what LoKI-B and BOLSIG+ assume for a bulk gas. Labels that are not
resolvable formulas still parse and round-trip — they just have no mass.

If [Catalyst.jl](https://catalyst.sciml.ai/stable/) is loaded, reactions can be
turned into Catalyst reactions that can be added to a modeling toolkit:

```julia
julia> using Catalyst

julia> to_catalyst(tree, (1.0, reaction))
1-element Vector{Reaction}:  
 Reaction(1.0, t, [e(t), N₂(t)], [2e(t), N₂⁺(t)])
```

## Project status

The package is under active development.  If you encounter an issue or have a
feature request, please open an [issue on GitHub](https://github.com/jqfeld/PlasmaSpecies.jl/issues).
