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
Species(StringGas("e"), Negative(1), nothing, nothing, nothing)

julia> nitrogen_ion = Species("N2(+,B)")
Species(Nitrogen(), Positive(1), StringElectronicState("B"), nothing, nothing)

julia> gas(nitrogen_ion), charge(nitrogen_ion)
(Nitrogen(), Positive(1))

julia> reaction = p"e + N2 --> 2e + N2(+)"
PlasmaReaction(Species[Species(StringGas("e"), Negative(1), nothing, nothing, nothing), Species(Nitrogen(), Neutral(), nothing, nothing, nothing)], Species[Species(StringGas("e"), Negative(1), nothing, nothing, nothing), Species(Nitrogen(), Positive(1), nothing, nothing, nothing)], [1, 1], [2, 1], false)

julia> reaction.subs, reaction.prods
(Species[Species(StringGas("e"), Negative(1), nothing, nothing, nothing), Species(Nitrogen(), Neutral(), nothing, nothing, nothing)], Species[Species(StringGas("e"), Negative(1), nothing, nothing, nothing), Species(Nitrogen(), Positive(1), nothing, nothing, nothing)])

julia> tree = SpeciesTree([electron, nitrogen_ion])
SpeciesTree
└─ Species(StringGas("e"), Negative(1), nothing, nothing, nothing)
└─ Species(Nitrogen(), Positive(1), StringElectronicState("B"), nothing, nothing)
```

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
