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

# Plain ASCII, which LoKI-B string round-tripping relies on. Isotope segments
# are '-'-joined per the ExoMol iso_slug spelling ("16O-1H"); "16O1H" would be
# ambiguous.
function Base.print(io::IO, g::Molecule)
  isotopic = any(p -> first(p) isa Isotope, g.composition)
  first_seg = true
  for (nuc, n) in g.composition
    isotopic && !first_seg && print(io, '-')
    print(io, nuc)
    n == 1 || print(io, n)
    first_seg = false
  end
end

function Base.show(io::IO, g::Molecule)
  for (nuc, n) in g.composition
    show(io, nuc)
    n == 1 || print(io, subscript(n))
  end
end

# A formula segment: optional mass number, element symbol, optional count.
const SEGMENT_RE = r"^(\d+)?([A-Z][a-z]?)(\d+)?$"
# Repeated-segment scan for a plain (non-isotopic) formula.
const PLAIN_SEGMENT_RE = r"([A-Z][a-z]?)(\d+)?"

"""
    parse_formula(s) -> Union{Molecule,Nothing}

Parse a chemical formula into a [`Molecule`](@ref), or return `nothing` if it
does not resolve — unknown element symbol, unknown isotope, or leftover text.
Returning `nothing` rather than throwing lets [`Gas`](@ref) fall back to
`StringGas`.

Three spellings are accepted:

- **Natural abundance**: `"Ar"`, `"N2"`, `"CO2"`, `"H2O"`, `"CH4"`. Every
  symbol resolves to an [`Element`](@ref), i.e. the standard atomic weight.
- **A single isotopic segment**: `"14N2"`, `"16O"`.
- **ExoMol `iso_slug`**: `"16O-1H"`, `"12C-16O2"`, `"28Si-16O"` — segments
  joined by `-`. Matches the spelling `ExoMol.jl` uses, so species strings line
  up with `load_isotopologue` slugs.

Two-letter symbols match greedily, per the usual chemical convention: `"CO"` is
carbon monoxide while `"Co"` is cobalt, `"NO"` is N + O while `"No"` is
nobelium. Case is significant.
"""
function parse_formula(s::AbstractString)
  isempty(s) && return nothing

  segments = if occursin('-', s)
    split(s, '-')
  elseif isdigit(first(s))
    [s]
  else
    return parse_plain_formula(s)
  end

  nucs = Nuclide[]
  counts = Int[]
  for seg in segments
    m = match(SEGMENT_RE, seg)
    isnothing(m) && return nothing
    sym = Symbol(m[2])
    nuc = if isnothing(m[1])
      element(sym)
    else
      isotope(sym, parse(Int, m[1]))
    end
    isnothing(nuc) && return nothing
    push!(nucs, nuc)
    push!(counts, isnothing(m[3]) ? 1 : parse(Int, m[3]))
  end
  return Molecule(_pairs_tuple(nucs, counts))
end

# `"CO2"` -> C, O2. Scanned rather than split, so the whole string has to be
# consumed; anything left over means the label is not a formula.
function parse_plain_formula(s::AbstractString)
  nucs = Nuclide[]
  counts = Int[]
  pos = 1
  while pos <= lastindex(s)
    m = match(PLAIN_SEGMENT_RE, s, pos)
    (isnothing(m) || m.offset != pos) && return nothing
    nuc = element(Symbol(m[1]))
    isnothing(nuc) && return nothing
    push!(nucs, nuc)
    push!(counts, isnothing(m[2]) ? 1 : parse(Int, m[2]))
    pos = m.offset + lastindex(m.match)
  end
  isempty(nucs) && return nothing
  return Molecule(_pairs_tuple(nucs, counts))
end
