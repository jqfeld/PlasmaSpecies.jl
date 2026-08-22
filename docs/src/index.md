```@meta
CurrentModule = PlasmaSpecies
```

# PlasmaSpecies.jl

Bookkeeping for the species and reactions of a plasma-chemistry model.

[PlasmaSpecies.jl](https://github.com/jqfeld/PlasmaSpecies.jl) gives one
concrete answer to "what is a species, and when are two of them the same" —
then builds the hierarchy, reaction parsing and export on top of it. It does no
physics: no rate coefficients, no transport, no solver.

## What it provides

- **[`Species`](@ref)** — a gas, a charge, and optionally an electronic,
  vibrational and rotational state. Written in the string notation of the
  [LoKI-B](https://github.com/IST-Lisbon/LoKI) Boltzmann solver
  (`"N2[+,B,vib=3]"`), with a `p"..."` macro for both species and reactions.
- **Gases resolved from chemical formulas** — anything built from known
  nuclides works without registration and knows its own mass, natural-abundance
  (`"CO2"`) or isotopically specific (`"12C-16O2"`).
- **[`SpeciesTree`](@ref)** — the parent/child hierarchy relating bulk `N2` to
  `N2[X]` to `N2[X,vib=1]`, with the leaves as the states a model integrates.
- **[`ReactionFormula`](@ref) / [`PlasmaReaction`](@ref)** — parsed reactions,
  reversible and fractional coefficients included, and
  [`apply_tree`](@ref) to expand one written over lumped species into one per
  resolved state.
- **[`to_catalyst`](@ref)** — export to
  [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/) reaction networks, via a
  package extension.

## Installation

```julia
julia> import Pkg; Pkg.add(url="https://github.com/jqfeld/PlasmaSpecies.jl.git")
```

## A first look

```@repl index
using PlasmaSpecies
sp = p"N2[+,B]"
gas(sp), charge(sp), electronic_state(sp)
mass(sp)
r = p"e + N2 --> 2e + N2[+]"
t = SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[+]"])
leaves(t)
apply_tree(t, p"e + N2[X] --> 2e + N2[+]")
```

## Where to go next

| Page | Covers |
|---|---|
| [Species](species.md) | the string notation, fields, identity, state matching |
| [Gases and masses](gases.md) | formula parsing, isotopes, nuclide data |
| [Species trees](trees.md) | building and querying the hierarchy |
| [Reactions](reactions.md) | parsing, expansion over a tree, energetics, Catalyst |
| [API](api.md) | every exported name |
