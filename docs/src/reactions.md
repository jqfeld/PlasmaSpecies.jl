```@meta
CurrentModule = PlasmaSpecies
```

# Reactions

Two types. A [`ReactionFormula`](@ref) is stoichiometry only — which species,
which coefficients, reversible or not. A [`PlasmaReaction`](@ref) is a formula
with a rate attached.

## Parsing

```@repl reactions
using PlasmaSpecies
r = p"e + N2 --> 2e + N2[+]"
r.subs, r.substoich
r.prods, r.prodstoich
r.reverse
```

The `p"..."` macro returns a `ReactionFormula` when the string contains an
arrow and a [`Species`](@ref) otherwise; [`parse_reaction`](@ref) is the
function form. Terms are separated by `+`, and a `+` inside brackets (`N2[+]`)
is not a separator.

Three arrows:

```@repl reactions
p"N2[X,vib=0] <--> N2[X,vib=1]"     # reversible, sets r.reverse
parse_reaction("N2[X,vib=1] <-- N2[X,vib=0]")   # normalised to forward
```

Coefficients may be fractional, which is common for wall and surface
reactions. A decimal anywhere on a side promotes that side's coefficients to
`Float64`; otherwise they stay `Int`:

```@repl reactions
f = parse_reaction("O[X] + 0.5O2[X] --> O3")
f.substoich
```

Both sides are normalised on construction: repeated species are combined and
species are sorted, so formulas that mean the same thing print and compare the
same.

```@repl reactions
parse_reaction("e + e + N2 --> N2")
```

## Expanding over a tree

A chemistry is usually written over lumped species and solved over resolved
ones. [`apply_tree`](@ref) does the expansion: every species in the formula is
resolved to the leaves it matches, and one reaction is emitted per combination.

```@repl reactions
t = SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]", "N2[+]"])
apply_tree(t, p"e + N2[X,vib=0-1] --> e + N2[X,vib=2]")
```

Each result carries a branching factor of `1 / (number of product
combinations)`: the rate is split equally over the possible outcomes. Writing a
product as a lumped species therefore distributes it:

```@repl reactions
apply_tree(t, p"e + N2[X,vib=0] --> e + N2[X]")
```

For a [`PlasmaReaction`](@ref) the factor is folded into the rate by
[`scale_rate`](@ref) instead, and the result is a vector of `PlasmaReaction`:

```@repl reactions
pr = PlasmaReaction(1.5e-16, p"e + N2[X,vib=0-1] --> e + N2[X,vib=2]")
[(x.rate, x.formula) for x in apply_tree(t, pr)]
```

A species that matches nothing in the tree would silently expand the reaction
away, so instead the reaction is dropped with a warning. Passing several
reactions concatenates the results and skips any that throw.

Expansion can produce identities (`A --> A`) when a lumped reaction collapses
onto one state. [`isreactive`](@ref) filters them out.

## Rates

`PlasmaReaction.rate` is deliberately untyped. A number is a constant rate
coefficient; anything callable is left to the consumer to evaluate, so a fit in
gas temperature or reduced field, or an interpolation from a Boltzmann solver,
works unchanged. [`scale_rate`](@ref) is the only place this package touches a
rate, and it handles both:

```@repl reactions
scale_rate(2.0, 0.5)
scale_rate(Te -> 1e-15 * Te^0.5, 0.5)(4.0)
```

## Energetics

Give the species energies — usually with [`apply_energy!`](@ref) over a tree —
and a reaction's energy balance follows:

```@repl reactions
tv = SpeciesTree(["N2[X,vib=0]", "N2[X,vib=1]"]);
apply_energy!(tv, q -> 0.29 * q.vib);
lv = leaves(tv)
pr = PlasmaReaction(2.0, ReactionFormula([lv[1]], [lv[2]], [1], [1], false))
reaction_energy(pr)
thermal_source_term(pr)(1e20)
```

[`reaction_energy`](@ref) is `Σ E(products) - Σ E(substrates)`, negative for an
exothermic reaction, and `nothing` if any energy is unset.
[`thermal_source_term`](@ref) builds the gas-heating term
`rate * ∏ nᵢ^stoichᵢ * (-ΔE)` as a function of the substrate densities, in the
order of `pr.formula.subs`.

Two consistency checks are available, though not exported:
`PlasmaSpecies.ismassbalanced` and `PlasmaSpecies.ischargebalanced`.

```@repl reactions
PlasmaSpecies.ischargebalanced(p"e + N2 --> 2e + N2[+]")
```

## Export to Catalyst

Loading [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/) activates an
extension providing [`to_catalyst`](@ref). Species variables are named after the
species string, so anything with a charge or state label becomes a `var"..."`
identifier.

```julia
julia> using PlasmaSpecies, Catalyst

julia> t = SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[+]"]);

julia> rs = to_catalyst(t, (1.0, p"e + N2[X] --> 2e + N2[+]"))
2-element Vector{Reaction{...}}:
 1.0, e + var"N2[X,vib=0]" --> 2*e + var"N2[+]"
 1.0, e + var"N2[X,vib=1]" --> 2*e + var"N2[+]"

julia> ReactionSystem(rs, default_t(); name = :plasma)
Model plasma:
Unknowns (4): see unknowns(plasma)
  e(t)
  var"N2[X,vib=0]"(t)
  var"N2[+]"(t)
  var"N2[X,vib=1]"(t)
```

Reactions are given as `(rate, formula)` tuples or `PlasmaReaction`s, singly, as
varargs or as a vector, and are expanded over the tree on the way out. A
`SpeciesTree` on its own converts to the sum of its leaf variables.

```@docs; canonical=false
ReactionFormula
PlasmaReaction
parse_reaction
@p_str
apply_tree
scale_rate
isreactive
reaction_energy
thermal_source_term
to_catalyst
```
