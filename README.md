# PlasmaSpecies

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jqfeld.github.io/PlasmaSpecies.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/PlasmaSpecies.jl/dev/)
[![Build Status](https://github.com/jqfeld/PlasmaSpecies.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jqfeld/PlasmaSpecies.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jqfeld/PlasmaSpecies.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jqfeld/PlasmaSpecies.jl)

Species and reaction bookkeeping for plasma-chemistry models.

`PlasmaSpecies.jl` defines what a species *is* — a gas, a charge, and
optionally an electronic, vibrational and rotational state — and builds the
hierarchy, reaction parsing and export on top of that. It does no physics:
no rate coefficients, no transport, no solver.

- Species and reactions in the string notation of
  [LoKI-B](https://github.com/IST-Lisbon/LoKI), with a `p"..."` macro.
- Gases resolved from chemical formulas, natural-abundance or isotopically
  specific, each knowing its own mass.
- A `SpeciesTree` relating bulk `N2` to `N2[X]` to `N2[X,vib=1]`, and
  `apply_tree` to expand a lumped reaction into one per resolved state.
- Export to [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/) reaction
  networks through a package extension.

Full documentation: <https://jqfeld.github.io/PlasmaSpecies.jl/stable/>

## Installation

```julia
julia> import Pkg; Pkg.add(url="https://github.com/jqfeld/PlasmaSpecies.jl.git")
```

## Quick start

```julia
julia> using PlasmaSpecies

julia> electron = p"e"
e

julia> nitrogen_ion = Species("N2[+,B]")
N2[+,B]

julia> gas(nitrogen_ion), charge(nitrogen_ion), electronic_state(nitrogen_ion)
(N₂, +, B)

julia> reaction = p"e + N2 --> 2e + N2[+]"
e+N2-->2e+N2[+]

julia> string.(reaction.prods), reaction.prodstoich
(["e", "N2[+]"], [2, 1])

julia> tree = SpeciesTree([electron, nitrogen_ion])
nothing
├─ e
└─ N2[+]
   └─ N2[+,B]
```

Species are compared and hashed on identity alone — gas, charge, and the three
state fields — so attaching an `energy`, `degeneracy` or `metadata` payload
never breaks a `Dict{Species,V}` lookup.

### Gases are compositions

A gas is resolved from its chemical formula, so anything built from known
nuclides works without being registered first, and knows its own mass:

```julia
julia> mass(Species("CO2")) / 1.66053906660e-27   # daltons
44.009

julia> gas(Species("H2O")), gas(Species("SF6"))
(H₂O, SF₆)
```

A bare element symbol means natural isotopic abundance (`N2` is 2 × 14.007 u),
which is what LoKI-B and BOLSIG+ assume for a bulk gas. Write a specific
isotopologue with ExoMol's `iso_slug` spelling, which lines these strings up
with `ExoMol.jl`'s `load_isotopologue`:

```julia
julia> mass(Species("28Si-16O")) / 1.66053906660e-27
43.97184115422001

julia> gas(Species("16O-1H"))
¹⁶O¹H
```

Labels that are not resolvable formulas still parse and round-trip — they just
have no mass.

### Reactions over a tree

A reaction written for a lumped species expands to the states actually present:

```julia
julia> t = SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[+]"]);

julia> apply_tree(t, p"e + N2[X] --> 2e + N2[+]")
2-element Vector{Tuple{Float64, ReactionFormula}}:
 (1.0, e+N2[X,vib=0]-->2e+N2[+])
 (1.0, e+N2[X,vib=1]-->2e+N2[+])
```

Each result carries a branching factor — one over the number of product
combinations — so a lumped *product* distributes the rate over its states.
Pair a formula with a rate as a `PlasmaReaction` and the factor is folded into
the rate instead.

### Catalyst export

Loading Catalyst activates the extension:

```julia
julia> using Catalyst

julia> to_catalyst(t, (1.0, reaction))
2-element Vector{Reaction{...}}:
 1.0, e + var"N2[X,vib=0]" --> 2*e + var"N2[+]"
 1.0, e + var"N2[X,vib=1]" --> 2*e + var"N2[+]"
```

Catalyst variables are named after the species string, so anything carrying a
charge or a state label becomes a `var"..."` identifier (`var"N2[+]"`,
`var"N2[X,vib=1]"`).

## Project status

Under active development. Issues and feature requests:
<https://github.com/jqfeld/PlasmaSpecies.jl/issues>.
