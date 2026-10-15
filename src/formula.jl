"""
    Molecule(composition) <: Gas

A gas defined by what it is made of: a tuple of `nuclide => count` pairs in the
order they were written. Mass and both display forms derive from the
composition, so any formula that resolves against the nuclide tables is a
usable gas.

`print` preserves the written order (`H2O`, not `OH2`), while `==` and `hash`
use the canonical composition, so `OH` and `HO` are one gas and one `Dict` key.
"""
struct Molecule{C<:Tuple,K<:Tuple} <: Gas
  composition::C
  canonical::K
end

# Canonicalise once here: `Species` is used as a `Dict` key downstream, so `==`
# and `hash` on a gas are on a hot path and must not allocate.
Molecule(composition::Tuple) = Molecule(composition, _canonicalize(composition))

function _canonicalize(composition::Tuple)
  nucs = Nuclide[]
  counts = Int[]
  for (nuc, n) in composition
    i = findfirst(==(nuc), nucs)
    if isnothing(i)
      push!(nucs, nuc)
      push!(counts, n)
    else
      counts[i] += n
    end
  end
  order = sortperm(nucs)
  return _pairs_tuple(nucs[order], counts[order])
end

# Built from the concrete values so the tuple gets concrete `Pair{Element,Int}` /
# `Pair{Isotope,Int}` element types rather than an abstract `Pair{Nuclide,Int}`.
_pairs_tuple(nucs, counts) = ntuple(i -> nucs[i] => counts[i], length(nucs))

"""
    composition(g) -> Tuple of `nuclide => count` pairs

Written-order composition of a gas. Errors for gases that have no resolved
composition (`StringGas`).
"""
composition(g::Molecule) = g.composition

"""
    canonical_composition(g) -> Tuple of `nuclide => count` pairs

Composition with repeated nuclides merged and sorted by (atomic number, mass
number). This is the identity used by `==` and `hash`.
"""
canonical_composition(g::Molecule) = g.canonical

mass(g::Molecule) = sum(n * mass(nuc) for (nuc, n) in g.composition)

Base.:(==)(a::Molecule, b::Molecule) = a.canonical == b.canonical
Base.hash(g::Molecule, h::UInt) = hash(g.canonical, hash(:Molecule, h))

# Plain ASCII, which LoKI-B string round-tripping relies on. A nuclide prints
# its mass number in braces ("{16}O{1}H"), the same spelling `parse_formula`
# reads, so no separator is needed between segments.
function Base.print(io::IO, g::Molecule)
  for (nuc, n) in g.composition
    print(io, nuc)
    n == 1 || print(io, n)
  end
end

function Base.show(io::IO, g::Molecule)
  for (nuc, n) in g.composition
    show(io, nuc)
    n == 1 || print(io, subscript(n))
  end
end

# One formula segment: optional braced mass number, element symbol, optional
# count. Scanned repeatedly from the start of the string.
const SEGMENT_RE = r"(?:\{(\d+)\})?([A-Z][a-z]?)(\d+)?"
# The pre-0.3 isotope spelling (a bare mass number: "14N2", "16O-1H"), which
# `Gas` rejects with a pointer to the new one rather than silently returning a
# massless `StringGas`.
const LEGACY_ISOTOPE_RE = r"^(?:\d*[A-Z][a-z]?\d*)(?:-\d*[A-Z][a-z]?\d*)*$"

"""
    parse_formula(s) -> Union{Molecule,Nothing}

Parse a chemical formula into a [`Molecule`](@ref), or return `nothing` if it
does not resolve — unknown element symbol, unknown isotope, or leftover text.
Returning `nothing` rather than throwing lets [`Gas`](@ref) fall back to
`StringGas`.

A formula is a run of segments, each an element symbol with an optional mass
number in braces before it and an optional count after it:

- **Natural abundance**: `"Ar"`, `"N2"`, `"CO2"`, `"H2O"`, `"CH4"`. Every
  symbol resolves to an [`Element`](@ref), i.e. the standard atomic weight.
- **Specific nuclides**: `"{14}N2"`, `"{16}O{1}H"`, `"{12}C{16}O2"`. The
  braced number is the mass number of the symbol that follows.
- **Mixed**: `"{13}CO2"` — one labelled atom, the rest at natural abundance.

The braces keep a mass number distinct from a stoichiometric coefficient, so
`2{16}O{1}H` in a reaction is unambiguous. For ExoMol's `iso_slug` spelling
(`"16O-1H"`) see [`from_iso_slug`](@ref) and [`iso_slug`](@ref).

Two-letter symbols match greedily, per the usual chemical convention: `"CO"` is
carbon monoxide while `"Co"` is cobalt, `"NO"` is N + O while `"No"` is
nobelium. Case is significant.
"""
function parse_formula(s::AbstractString)
  nucs = Nuclide[]
  counts = Int[]
  pos = 1
  while pos <= lastindex(s)
    m = match(SEGMENT_RE, s, pos)
    (isnothing(m) || m.offset != pos) && return nothing
    sym = Symbol(m[2])
    nuc = isnothing(m[1]) ? element(sym) : isotope(sym, parse(Int, m[1]))
    isnothing(nuc) && return nothing
    push!(nucs, nuc)
    push!(counts, isnothing(m[3]) ? 1 : parse(Int, m[3]))
    pos = m.offset + lastindex(m.match)
  end
  isempty(nucs) && return nothing
  return Molecule(_pairs_tuple(nucs, counts))
end

# Whether `s` is spelled the pre-0.3 way: a bare mass number in front of a
# symbol, segments optionally '-'-joined. Only strings that would have resolved
# as isotopes under the old grammar qualify, so ordinary labels pass through.
is_legacy_isotope_spelling(s::AbstractString) =
  occursin(LEGACY_ISOTOPE_RE, s) && occursin(r"(?:^|-)\d+[A-Z]", s)

"""
    from_iso_slug(slug) -> Molecule

Build a gas from an ExoMol `iso_slug` — segments of mass number, symbol and
optional count joined by `-`, as `ExoMol.jl`'s `load_isotopologue` takes them:
`"16O-1H"`, `"12C-16O2"`, `"28Si-16O"`. Every segment must carry a mass
number. Throws `ArgumentError` for a malformed slug or an unknown nuclide.
Inverse of [`iso_slug`](@ref).
"""
function from_iso_slug(slug::AbstractString)
  nucs = Nuclide[]
  counts = Int[]
  for seg in split(slug, '-')
    m = match(r"^(\d+)([A-Z][a-z]?)(\d+)?$", seg)
    isnothing(m) && throw(ArgumentError(
      "$(repr(slug)) is not an ExoMol iso_slug: segment $(repr(seg)) should be <mass number><symbol>[<count>]"))
    nuc = isotope(Symbol(m[2]), parse(Int, m[1]))
    isnothing(nuc) && throw(ArgumentError("unknown nuclide $(repr(seg)) in iso_slug $(repr(slug))"))
    push!(nucs, nuc)
    push!(counts, isnothing(m[3]) ? 1 : parse(Int, m[3]))
  end
  return Molecule(_pairs_tuple(nucs, counts))
end

"""
    iso_slug(g::Molecule) -> String

The ExoMol `iso_slug` spelling of a fully isotope-resolved gas, in written
order: `{16}O{1}H` gives `"16O-1H"`. Throws `ArgumentError` if any nuclide is
at natural abundance, because a slug names one isotopologue. Inverse of
[`from_iso_slug`](@ref).
"""
function iso_slug(g::Molecule)
  any(p -> first(p) isa Element, g.composition) && throw(ArgumentError(
    "$(g) has natural-abundance nuclides; an iso_slug names a single isotopologue"))
  join(("$(nuc.A)$(nuc.symbol)$(n == 1 ? "" : n)" for (nuc, n) in g.composition), '-')
end
