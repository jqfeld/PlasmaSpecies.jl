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

[`parse_formula`](@ref) accepts three spellings.

**Natural abundance.** Bare element symbols, optionally with counts. Each
resolves to an [`Element`](@ref), i.e. the standard atomic weight — `N2` is
2 × 14.007 u, which is what LoKI-B and BOLSIG+ assume for a bulk gas.

```@repl gases
composition(gas(Species("CH4")))
```

**A single isotopic segment.** A mass number in front: `"14N2"`, `"16O"`.

**ExoMol `iso_slug`.** Segments joined by `-`, matching the spelling
`ExoMol.jl` uses, so species strings line up with `load_isotopologue` slugs.

```@repl gases
gas(Species("16O-1H"))
gas(Species("12C-16O2"))
string(gas(Species("16O-1H")))
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
species parser only accepts letters, digits and `-` in the gas token, so a
label like `"my_gas"` is a parse error, while `"mygas"` is a `StringGas`.

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
mass(Species("28Si-16O")) / PlasmaSpecies.U_KG
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
