```@meta
CurrentModule = PlasmaSpecies
```

# Gases and masses

The [`Gas`](@ref) of a [`Species`](@ref) is its chemical identity, independent
of charge and internal state. A gas label is resolved by `Gas(label)` in three
steps:

1. [`REGISTERED_GASES`](@ref) — non-formula labels. Only `"e"`, giving
   [`Electron`](@ref).
2. [`parse_formula`](@ref) — a chemical formula, giving a [`Molecule`](@ref)
   that knows its composition and therefore its mass.
3. [`StringGas`](@ref) — a label that resolves to neither. It parses, prints and
   round-trips, but has no composition and no mass.

There is no registration step: any formula built from the nuclide tables is a
usable gas.

```@repl gases
using PlasmaSpecies
gas(Species("CO2"))
gas(Species("SF6"))
gas(Species("Xyz"))
```

## Formula syntax

[`parse_formula`](@ref) reads a run of segments, each an element symbol with an
optional mass number in braces before it and an optional count after it.

**Natural abundance.** Bare element symbols, optionally with counts. Each
resolves to an [`Element`](@ref), i.e. the standard atomic weight — `N2` is
2 × 14.007 u, which is what LoKI-B and BOLSIG+ assume for a bulk gas.

```@repl gases
composition(gas(Species("CH4")))
```

**Specific nuclides.** A mass number in braces in front of the symbol it
labels: `"{14}N2"`, `"{16}O{1}H"`, `"{12}C{16}O2"`. Labelled and
natural-abundance atoms mix freely: `"{13}CO2"`.

```@repl gases
gas(Species("{16}O{1}H"))
gas(Species("{13}CO2"))
string(gas(Species("{16}O{1}H")))
```

The braces keep a mass number distinct from a stoichiometric coefficient, so
`2{16}O{1}H` in a reaction is two ¹⁶O¹H. The pre-0.3 spelling with a bare
mass number (`"14N2"`, `"16O-1H"`) is rejected with an error pointing here.

**ExoMol interop.** [`iso_slug`](@ref) and [`from_iso_slug`](@ref) convert to
and from ExoMol's `iso_slug` spelling (`"16O-1H"`), the form `ExoMol.jl`'s
`load_isotopologue` takes.

```@repl gases
iso_slug(gas(Species("{12}C{16}O2")))
from_iso_slug("28Si-16O")
```

Two display forms: `print`/`string` is plain ASCII and preserves the written
order, so a species round-trips through its string; `show` is the pretty form
with sub- and superscripts.

Two-letter symbols match greedily, per the usual chemical convention, and case
is significant:

```@repl gases
gas(Species("CO"))   # carbon monoxide
gas(Species("Co"))   # cobalt
```

A formula that does not resolve — unknown symbol, unknown isotope, leftover
text — falls back to [`StringGas`](@ref) rather than throwing. Note that the
species parser only accepts letters, digits, braces and `-` in the gas token,
so a label like `"my_gas"` is a parse error, while `"mygas"` is a `StringGas`.

## Composition and equality

A [`Molecule`](@ref) stores the composition twice: as written, and canonical
(repeated nuclides merged, sorted by atomic and mass number). `==` and `hash`
use the canonical form, so `OH` and `HO` are the same gas and the same `Dict`
key, while display keeps what was written.

```@repl gases
gas(Species("OH")) == gas(Species("HO"))
composition(gas(Species("H2O")))
canonical_composition(gas(Species("H2O")))
```

## Masses

[`mass`](@ref) is in kilograms throughout. For a species, the electrons implied
by the charge are added or removed.

```@repl gases
mass(Species("CO2")) / PlasmaSpecies.U_KG   # in daltons
mass(Species("{28}Si{16}O")) / PlasmaSpecies.U_KG
mass(Species("N2[+]")) == mass(Species("N2")) - mass(Species("e"))
```

## Nuclide data

Two tables, both in `src/nuclide_data.jl`:

- [`ELEMENTS`](@ref) — all 118 elements, standard atomic weights (IUPAC/CIAAW
  2021 conventional values). Elements with no stable nuclide carry the mass
  number of their longest-lived isotope.
- [`ISOTOPES`](@ref) — exact nuclide masses (AME2020), keyed by
  `(symbol, mass number)`. Not the full nuclide chart: it covers the elements
  that come up in plasma chemistry and molecular spectroscopy. Adding one is a
  single line.

[`ELEMENT_ALIASES`](@ref) maps symbols that stand for a specific nuclide —
`D` and `T` — so `"D2"` and `"HT"` resolve without a mass number.

```@repl gases
PlasmaSpecies.element(:D)
mass(Species("D2")) / PlasmaSpecies.U_KG
```

Look-ups are [`element`](@ref) and [`isotope`](@ref), both returning `nothing`
for anything unknown.

## Gases without a formula

A [`StringGas`](@ref) behaves like any other gas except that it cannot be
weighed:

```@repl gases
mass(Species("Xyz"))
```

If a label needs a mass anyway, define a gas type and a `mass` method for it:

```@repl gases
struct UnknownGas <: PlasmaSpecies.Gas end
PlasmaSpecies.mass(::UnknownGas) = π * 1e-26
mass(Species(UnknownGas()))
```
