# Changelog

## 0.3.0

### Breaking

- **Mass numbers are written in braces** before the symbol they label:
  `{14}N2`, `{16}O{1}H`, `{12}C{16}O2`, `{13}CO2`. The bare-digit and ExoMol
  slug spellings (`14N2`, `16O-1H`) are rejected by `Gas`/`Species` with an
  `ArgumentError` naming the new form. Leading digits in a reaction term are
  therefore always the stoichiometric coefficient (`2{16}O{1}H`).
- ExoMol interop moves to explicit helpers: `iso_slug(::Molecule)` and
  `from_iso_slug(str)`. Isotopes print as `{A}Sym`, so `D2` prints `{2}H2`.
- `ReactionFormula{S,P}` is parametrised on the coefficient type of each side
  (`ReactionFormula{Int,Float64}` for `O[3P] --> 0.5O2[X]`), with concrete
  `Vector{Species}` / `Vector{S}` / `Vector{P}` fields instead of untyped ones.
  Promotion is per side, as before. The constructor rejects a side whose
  species and coefficient counts differ.
- The batch forms of `apply_tree` and `to_catalyst` no longer catch errors.
  An error in any reaction is logged with that reaction and rethrown, aborting
  the batch. They return typed vectors
  (`Vector{Tuple{Float64,ReactionFormula}}`, `Vector{PlasmaReaction}`,
  `Vector{Reaction}`) instead of `Vector{Any}`.
- Species of one gas sort with resolved states before unresolved ones at every
  level, and the ground state `X` first: `N2[X,vib=0] < N2[X] < N2[A] < N2`.
- Julia 1.10 is the minimum; `[compat]` bounds added (Catalyst 15–16).

### Added

- `thermal_source_term` evaluates a callable rate: the returned function takes
  `rate_args`, forwarded as `rate(rate_args...)`, so one closure follows a
  changing electron temperature.
- `isnegative`/`ispositive` extend `Base`'s on Julia ≥ 1.13 instead of
  shadowing them.
- `SpeciesTree()` builds an empty tree; `Species` and `SpeciesTree` accept any
  `AbstractString`.
- `[0]` is an explicit neutral charge.

### Fixed

- Equal `ElectronicState`s with different representations now hash equal, so
  `Species("N2[B3Pig]")` and `Species("N2[B]")` are one `Dict` key.
- `isless` on species with equal identity fields returned `nothing`; sorting a
  vector with duplicates threw.
- Multiply charged ions (`Ar[+++]`) parse inside reactions.
- `N2[vib=1]` is a vibrational level without an electronic state, not an
  electronic state named `vib=1`; `rot=` likewise. Such species print in a form
  that parses back.
- 1-tuple quantum labels (`vib=(1,)`) round-trip through `string`.
- `parse_reaction` without an arrow throws an `ArgumentError` instead of a
  `MethodError`.
- Whitespace stripping no longer throws on non-ASCII input.
- The Catalyst extension uses `default_t()`, so exported species share the
  independent variable of any `ReactionSystem`.

### CI

- Tests run on Julia LTS, the latest release and nightly, plus a job with every
  dependency at its compat lower bound.
