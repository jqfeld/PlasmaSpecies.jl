```@meta
CurrentModule = PlasmaSpecies
```

# Species and states

A [`Species`](@ref) is one state of one gas: a [`Gas`](@ref), a
[`Charge`](@ref), and optionally an electronic, vibrational and rotational
state.

```@docs; canonical=false
Species
```

## The string notation

Species are normally written as strings, in the notation of the
[LoKI-B](https://github.com/IST-Lisbon/LoKI) Boltzmann solver:

    GAS[CHARGE,ELECTRONIC,vib=VIB,rot=ROT]

Every bracketed field is optional, but the order is fixed.

| Field | Syntax | Becomes |
|---|---|---|
| `GAS` | `N2`, `CO2`, `16O-1H`, `e` | a [`Gas`](@ref) — see [Gases and masses](gases.md) |
| `CHARGE` | `+`, `++`, `-`, `0`, omitted | [`Positive`](@ref), [`Negative`](@ref), [`Neutral`](@ref) |
| `ELECTRONIC` | `X`, `a''`, `B3Pig` | an [`ElectronicState`](@ref) |
| `vib=`, `rot=` | `0`, `0-4`, `(0,1,0)`, `(v1=1,v2=0)`, `asym` | a [`QuantumLabel`](@ref) |

```@repl species
using PlasmaSpecies
Species("e")
Species("N2[+]")
Species("N2[+,B]")
Species("N2[X,vib=0,rot=10]")
p"N4[+]"
```

The `p"..."` macro parses a species or a reaction, whichever the string looks
like. Whitespace is stripped, and a `+` inside brackets is never mistaken for a
charge separator.

A malformed string throws rather than parsing the part it recognises —
`"N2[+"` is a typo, not a neutral N₂.

```@repl species
Species("N2[+")
```

### Quantum labels

`vib=`/`rot=` accept more than a single quantum number. The forms are collected
in [`QuantumLabel`](@ref) and parsed by [`parse_quantum_label`](@ref):

```@repl species
vibrational_state(Species("N2[X,vib=3]"))          # one level
vibrational_state(Species("N2[X,vib=0-4]"))        # a lumped range
vibrational_state(Species("CO2[X,vib=(0,1,0)]"))   # normal modes
vibrational_state(Species("CO2[X,vib=(v1=1,v2=0)]"))
vibrational_state(Species("N2[X,vib=asym]"))       # opaque label
```

A range is the one form with containment semantics: `vib=0-4` matches the
levels inside it when a reaction is expanded over a tree. A `Symbol` label is
deliberately opaque and only ever matches itself.

Display round-trips: `Species(string(sp)) == sp` for every form above.

### Electronic states

A plain letter label becomes a [`StringElectronicState`](@ref). A fused
term-symbol string, in the spelling ExoMol's `ElecState` column uses, decomposes
into a [`TermSymbol`](@ref):

```@repl species
electronic_state(Species("N2[+,B]"))
ts = ElectronicState("A3Sigmau+")
ts.multiplicity, ts.Λ, ts.parity, ts.reflection
```

Both compare by [`label`](@ref) alone, so a `TermSymbol` and a bare label for
the same state are equal and interchangeable when matching. The extra physics is
annotation.

## Fields and accessors

```@repl species
sp = Species("N2[+,B]")
gas(sp)
charge(sp)
electronic_state(sp)
vibrational_state(sp) === nothing
```

`nothing` means *unresolved*, not *zero*. `N2` is bulk N₂ and `N2[X,vib=0]` is
its ground vibrational level; they are different species, related by a
[`SpeciesTree`](@ref).

Three further fields carry data *about* a state rather than naming it:
[`energy`](@ref), [`degeneracy`](@ref) and [`metadata`](@ref). Set them with
[`with_fields`](@ref), which copies a species and replaces only the named
fields:

```@repl species
excited = with_fields(sp, energy = 11.03, degeneracy = 3)
energy(excited), degeneracy(excited)
```

Species are immutable, so this is the only way to change one.

## Identity

`==` and `hash` use the identity fields only — gas, charge, and the three state
fields. `energy`, `degeneracy` and `metadata` are ignored:

```@repl species
sp == excited
hash(sp) == hash(excited)
d = Dict(sp => "N2+ in B"); d[excited]
```

That is deliberate. Attaching an energy or a downstream payload must not break
a `Dict{Species,V}` lookup, which is how consumers index their state vectors.

Species also sort (`isless` is defined for them): electrons last, then by charge, gas,
electronic state (with `X` first), and level. Sorting is total and never throws,
even for mixed label types.

## Matching states

Matching is separate from equality: a species with `nothing` or a range in a
field acts as a pattern for the more resolved states under it.

```@docs; canonical=false
level_matches
species_matches
Base.in
is_parent_species
```

```@repl species
Species("N2[X,vib=1]") in Species("N2[X]")        # wildcard field
Species("N2[X,vib=2]") in Species("N2[X,vib=0-4]") # range contains level
Species("N2[X]") in Species("N2[X]")               # reflexive
is_parent_species(Species("N2[X]"), Species("N2[X]"))  # strict: not its own parent
```

Use `in`/[`species_matches`](@ref) to ask "does this pattern cover that state",
and [`is_parent_species`](@ref) for hierarchy construction, where a species must
not be its own ancestor. [`get_parent_species`](@ref) walks one step up:

```@repl species
get_parent_species(Species("N2[X,vib=1]"))
get_parent_species(Species("N2"))
```

## Mass

[`mass`](@ref) returns kg, derived from the gas composition, with the electrons
implied by the charge added or removed:

```@repl species
mass(p"N2")
mass(p"N2[+]") == mass(p"N2") - mass(p"e")
mass(p"e")
```

Any formula built from the nuclide tables works — no registration step. A label
that is not a resolvable formula still parses and round-trips, but has no mass:

```@repl species
gas(Species("Xyz"))
mass(Species("Xyz"))
```

See [Gases and masses](gases.md) for formula syntax, isotopes and the data
tables.
