"""
    Species(; gas, charge=Neutral(), electronic_state=nothing,
              vibrational_state=nothing, rotational_state=nothing,
              energy=nothing, degeneracy=nothing, metadata=nothing)

One state of one gas. All fields are type parameters, so a `Species` is
concretely typed and cheap to store in arrays and `Dict` keys.

Usually built from a LoKI-B string instead — see [`Species(::String)`](@ref).

## Identity fields

Compared by `==` and hashed. Together they name the state:

- `gas`: the parent gas. Normally a [`Molecule`](@ref) resolved from a chemical
  formula (see [`parse_formula`](@ref)), otherwise [`Electron`](@ref) or the
  [`StringGas`](@ref) fallback.
- `charge`: [`Neutral`](@ref) (default), `Positive(n)` or `Negative(n)`.
- `electronic_state`: an [`ElectronicState`](@ref), or `nothing` for unresolved.
- `vibrational_state`, `rotational_state`: a [`QuantumLabel`](@ref), or
  `nothing`. By convention a resolved level implies its parent is resolved too
  (`rot` without `vib` is meaningless), but this is not enforced.

`nothing` means *unresolved*, not *zero*: `N2` is bulk N₂, `N2[X,vib=0]` is its
ground vibrational level, and the two are different species that a
[`SpeciesTree`](@ref) relates as parent and child.

## Payload fields

Ignored by `==` and `hash`, carried through [`with_fields`](@ref):

- `energy`: level energy, used by [`reaction_energy`](@ref). Set in bulk with
  [`apply_energy!`](@ref).
- `degeneracy`: statistical weight of the level.
- `metadata`: opaque payload for downstream packages.

The split is deliberate: `energy`/`degeneracy`/`metadata` are annotations on a
state, not part of its identity, so attaching them never breaks a
`Dict{Species,V}` lookup.

An [`Electron`](@ref) gas always carries `Negative(1)`; any charge passed for one
is normalised to that.
"""
struct Species{S,C<:Charge,E,V,J,En,D,MD}
  gas::S
  charge::C
  electronic_state::E
  vibrational_state::V
  rotational_state::J
  energy::En
  degeneracy::D
  metadata::MD
end

function Species(;
  gas,
  charge = Neutral(),
  electronic_state = nothing,
  vibrational_state = nothing,
  rotational_state = nothing,
  energy = nothing,
  degeneracy = nothing,
  metadata = nothing,
)
  gas isa Electron && (charge = Negative(1))
  return Species(gas, charge, electronic_state, vibrational_state,
                 rotational_state, energy, degeneracy, metadata)
end

"""
    with_fields(sp::Species; kwargs...) -> Species

Copy `sp`, replacing only the fields named in `kwargs`. Everything else —
`energy`, `degeneracy` and the opaque `metadata` payload downstream packages
attach — is carried over.
"""
with_fields(sp::Species; kwargs...) = Species(;
  NamedTuple{fieldnames(Species)}(ntuple(i -> getfield(sp, i), nfields(sp)))...,
  kwargs...,
)

"""
    QuantumLabel

The closed set of forms `vibrational_state`/`rotational_state` can take: a single
quantum number (`Int`), several quantum numbers (`Tuple`/`NamedTuple`, e.g. normal
modes of a polyatomic), a contiguous range of levels lumped together for reaction
purposes (`UnitRange{Int}`, see [`level_matches`](@ref)/[`matching_leaves`](@ref)),
or an opaque hand-picked label (`Symbol`) with no containment semantics.
"""
const QuantumLabel = Union{Int,NTuple{N,Int} where N,NamedTuple,UnitRange{Int},Symbol}

"""
    parse_quantum_label(s) -> QuantumLabel

Parse the string form of a `vibrational_state`/`rotational_state` field into a
[`QuantumLabel`](@ref). `"3"` becomes an `Int`, `"0-4"` a `UnitRange{Int}`,
`"(1,0,2)"` a `Tuple`, `"(v1=1,v2=0)"` a `NamedTuple`, and anything else an
opaque `Symbol`. `nothing` passes through unchanged. Inverse of
[`format_quantum_label`](@ref).
"""
parse_quantum_label(::Nothing) = nothing
function parse_quantum_label(s::AbstractString)
  i = tryparse(Int, s)
  isnothing(i) || return i

  m = match(r"^(\d+)-(\d+)$", s)
  isnothing(m) || return UnitRange(parse(Int, m[1]), parse(Int, m[2]))

  if startswith(s, "(") && endswith(s, ")")
    parts = strip.(split(s[2:end-1], ","))
    if !isempty(parts) && all(p -> occursin("=", p), parts)
      pairs = map(parts) do p
        k, v = split(p, "=")
        Symbol(strip(k)) => parse(Int, strip(v))
      end
      return (; pairs...)
    end
    ints = tryparse.(Int, parts)
    if !isempty(ints) && all(!isnothing, ints)
      return Tuple(ints)
    end
  end

  return Symbol(s)
end


"""
    isnegative(x) -> Bool

Whether `x` is a negatively charged `Species`. Dispatches on the charge type
parameter, so it resolves at compile time. `false` for a non-`Species`.
"""
isnegative(x) = false
isnegative(::Species{S,Negative,E,V,J}) where {S,E,V,J} = true

"""
    ispositive(x) -> Bool

Whether `x` is a positively charged `Species`. See [`isnegative`](@ref).
"""
ispositive(x) = false
ispositive(::Species{S,Positive,E,V,J}) where {S,E,V,J} = true

"""
    isneutral(x) -> Bool

Whether `x` is an uncharged `Species`. See [`isnegative`](@ref).
"""
isneutral(x) = false
isneutral(::Species{S,Neutral,E,V,J}) where {S,E,V,J} = true

Base.:(==)(x::S1, y::S2) where {S1<:Species, S2<:Species} =
  gas(x) == gas(y) &&
  charge(x) == charge(y) &&
  electronic_state(x) == electronic_state(y) &&
  vibrational_state(x) == vibrational_state(y) &&
  rotational_state(x) == rotational_state(y)

# Must mirror `==` exactly (same field set, deliberately excluding
# `energy`/`degeneracy`/`metadata`): Julia's `Dict`/`Set` require
# `a == b` to imply `hash(a) == hash(b)`. Without this, the default
# struct-derived `hash` (which *does* include every field) disagrees with
# the custom `==` above the moment two otherwise-identical species differ
# only in `metadata` (or `energy`/`degeneracy`) -- confirmed directly to
# silently break `Dict{Species,V}` lookups (`d[s1]=...; d[s2]` missing
# even though `s1 == s2`), exactly the mechanism
# `PlasmaFluidSim.compile_reactions`'s `species_index` Dict depends on
# staying correct. `0x5ea3c1e5` is just a fixed salt distinguishing this
# hash from an unrelated type that happened to hash its own fields the
# same way; any constant works.
const _SPECIES_HASH_SALT = UInt === UInt64 ? 0x5ea3c1e5b7d4a291 : 0x5ea3c1e5
function Base.hash(x::Species, h::UInt)
  h = hash(gas(x), h)
  h = hash(charge(x), h)
  h = hash(electronic_state(x), h)
  h = hash(vibrational_state(x), h)
  h = hash(rotational_state(x), h)
  return hash(_SPECIES_HASH_SALT, h)
end


"""
    Species(str::String; metadata=nothing) -> Species

Parse a species string in the notation of the
[LoKI-B](https://github.com/IST-Lisbon/LoKI) Boltzmann solver:

    GAS[CHARGE,ELECTRONIC,vib=VIB,rot=ROT]

Every bracketed field is optional, but the order is fixed. `GAS` is resolved by
[`Gas`](@ref) — a chemical formula, `"e"`, or an unresolvable label kept as a
[`StringGas`](@ref). `CHARGE` is a run of `+`/`-` (see [`Charge`](@ref)), the
electronic state is parsed by [`ElectronicState`](@ref), and the `vib=`/`rot=`
levels by [`parse_quantum_label`](@ref). Whitespace is stripped first.

```julia
Species("e")
Species("N2")                 # bulk N₂, no state resolved
Species("N2[+]")
Species("N2[+,B]")
Species("N2[X,vib=0,rot=10]")
Species("N2[X,vib=0-4]")      # a lumped range of levels
Species("CO2[X,vib=(0,1,0)]") # normal modes of a polyatomic
Species("16O-1H")             # a specific isotopologue
```

The gas token accepts letters, digits and `-`, so ExoMol-style isotope slugs
stay one token. Anything else — an underscore, a stray bracket — throws:
`"N2[+"` is a typo, not a neutral N₂.

`metadata` is stored as-is and never inspected.
"""
function Species(str::String; metadata=nothing)
  str = join(map(x -> isspace(str[x]) ? "" : str[x], 1:length(str)))
  # The gas group accepts '-' so ExoMol-style isotope slugs ("16O-1H",
  # "12C-16O2") stay a single token instead of being cut at the first dash.
  # Anchored at both ends so trailing or malformed text is an error.
  regex = r"(^[A-Za-z0-9\-]*)(?:\[([+\-]+)?,?([^,\]]+)?(?:,vib=(\([^\)]*\)|[^,\]]*)(?:,rot=(\([^\)]*\)|[^,\]]*))?)?\])?$"
  matched = match(regex, str)
  isnothing(matched) && error(
    """Cannot parse species string $(repr(str)). Expected \
    GAS[(CHARGE,ELECTRONIC_STATE,vib=VSTATE,rot=JSTATE)], e.g. "N2[+,B,vib=3,rot=1]".""")
  m = collect(matched)
  Species(;
    gas = Gas(m[1]),
    charge = Charge(m[2]),
    electronic_state = ElectronicState(m[3]),
    vibrational_state = parse_quantum_label(m[4]),
    rotational_state = parse_quantum_label(m[5]),
    metadata)
end
Species(gas::G) where {G<:Gas} = Species(gas=gas)


"""
    gas(sp::Species) -> Gas

The parent gas. See [`Gas`](@ref).
"""
gas(sp::Species) = sp.gas

"""
    charge(sp::Species) -> Charge

The charge state. See [`Charge`](@ref).
"""
charge(sp::Species) = sp.charge

"""
    electronic_state(sp::Species) -> Union{ElectronicState,Nothing}

The electronic state, or `nothing` if unresolved.
"""
electronic_state(sp::Species) = sp.electronic_state

"""
    vibrational_state(sp::Species) -> Union{QuantumLabel,Nothing}

The vibrational level, or `nothing` if unresolved.
"""
vibrational_state(sp::Species) = sp.vibrational_state

"""
    rotational_state(sp::Species) -> Union{QuantumLabel,Nothing}

The rotational level, or `nothing` if unresolved.
"""
rotational_state(sp::Species) = sp.rotational_state

"""
    energy(sp::Species)

Level energy, or `nothing` if unset. Not part of species identity; set it with
[`with_fields`](@ref) or, for a whole tree, [`apply_energy!`](@ref). Units are
whatever the caller used — [`reaction_energy`](@ref) only subtracts them.
"""
energy(sp::Species) = sp.energy

"""
    degeneracy(sp::Species)

Statistical weight of the level, or `nothing` if unset. Not part of species
identity.
"""
degeneracy(sp::Species) = sp.degeneracy

"""
    metadata(sp::Species)

Opaque payload attached by a downstream package, or `nothing`. Not part of
species identity, and never inspected by this package.
"""
metadata(sp::Species) = sp.metadata

"""
    mass(sp::Species) -> Float64

Mass in kg, derived from the gas composition with the charge accounted for:
`n` missing electrons for `Positive(n)`, `n` extra for `Negative(n)`. An
electron species returns the electron rest mass directly.

Errors if the gas is a [`StringGas`](@ref), i.e. if the label never resolved to
a chemical formula.
"""
mass(sp::Species) = gas(sp) isa Electron ? mass(gas(sp)) : mass(gas(sp)) - to_value(charge(sp)) * mass(Electron())

"""
    get_parent_species(sp::Species) -> Union{Species,Nothing}

Strip the most specific resolved state field and return the resulting species:
`rot` first, then `vib`, then the electronic state. Returns `nothing` for a
species that has none of them, i.e. the root of its branch.

Applied repeatedly by [`SpeciesTree`](@ref) to build the chain of ancestors
above a species. Payload fields (`energy`, `degeneracy`, `metadata`) are not
carried over — a parent is a different state.
"""
function get_parent_species(sp::Species)
  if !isnothing(rotational_state(sp))
    return Species(;
      gas = gas(sp),
      charge = charge(sp),
      electronic_state = electronic_state(sp),
      vibrational_state = vibrational_state(sp),
    )
  elseif !isnothing(vibrational_state(sp))
    return Species(
      gas = gas(sp),
      charge = charge(sp),
      electronic_state = electronic_state(sp),
    )
  elseif !isnothing(electronic_state(sp))
    return Species(
      gas = gas(sp),
      charge = charge(sp),
    )
  else
    return nothing
  end
end

"""
    level_matches(formula_value, candidate_value) -> Bool

Whether a single field (electronic/vibrational/rotational state) of a formula
`Species` matches a concrete candidate value. `nothing` on the formula side is a
wildcard (matches anything, as used by [`is_parent_species`](@ref) for hierarchy
matching); a `UnitRange{Int}` on the formula side matches any `Int` it contains
(used by [`matching_leaves`](@ref) to expand reactions written over a range of
vibrational levels); everything else falls back to exact equality — including
`Symbol` labels, which are deliberately opaque (no containment).
"""
level_matches(::Nothing, _) = true
level_matches(r::UnitRange{Int}, x::Int) = x in r
level_matches(a, b) = a == b

"""
    species_matches(formula::Species, candidate::Species) -> Bool

Whether `candidate` is matched by `formula`: `gas`/`charge` compared exactly,
electronic/vibrational/rotational state compared via [`level_matches`](@ref)
(so `nothing`/range fields on `formula` act as wildcards/containment checks).
Used by [`matching_leaves`](@ref) to expand a reaction formula species — which
may specify a range of vibrational levels — into every matching leaf of a
`SpeciesTree`.
"""
species_matches(formula::Species, candidate::Species) =
  gas(formula) == gas(candidate) &&
  charge(formula) == charge(candidate) &&
  level_matches(electronic_state(formula), electronic_state(candidate)) &&
  level_matches(vibrational_state(formula), vibrational_state(candidate)) &&
  level_matches(rotational_state(formula), rotational_state(candidate))

"""
    candidate::Species in formula::Species -> Bool

Whether `candidate` is matched by `formula` — an alias for
[`species_matches`](@ref)`(formula, candidate)`, exposed via `∈`/`in` to mirror
`SpectraUtils`'s `Base.in(child::State, parent::State)` (e.g. `s in
boltzmann.parent`). `nothing`/range fields on `formula` act as wildcards/
containment checks, same as `species_matches`; reflexive (`sp in sp` is
`true`), unlike [`is_parent_species`](@ref) which is a strict-ancestor check
and doesn't compare `rotational_state`.
"""
Base.in(candidate::Species, formula::Species) = species_matches(formula, candidate)

"""
    is_parent_species(sp::Species, parent::Species) -> Bool

Whether `parent` is a strict ancestor of `sp` in the species hierarchy: same
`gas` and `charge`, electronic and vibrational state matched via
[`level_matches`](@ref) (so `nothing` on `parent` acts as a wildcard), and
`parent != sp`. `rotational_state` is deliberately not compared. Unlike
[`species_matches`](@ref), this is irreflexive — a species is not its own parent.
"""
function is_parent_species(sp::Species, parent::Species)
  if parent != sp && gas(parent) == gas(sp) && charge(parent) == charge(sp)
    if level_matches(electronic_state(parent), electronic_state(sp))
      return level_matches(vibrational_state(parent), vibrational_state(sp))
    else
      false
    end
  else
    return false
  end
end


"""
    format_quantum_label(x) -> String

Inverse of [`parse_quantum_label`](@ref), used by `Base.show(::Species)` so
`vib=`/`rot=` round-trip through display back to the same parsed value. Only
needed for `UnitRange`, whose default `string` uses `"a:b"` (`0:4`) rather than
the `"a-b"` (`0-4`) syntax `parse_quantum_label` accepts; every other
`QuantumLabel` variant's default `string` is already what the parser expects.
"""
format_quantum_label(r::UnitRange{Int}) = "$(first(r))-$(last(r))"
format_quantum_label(x) = string(x)

function Base.show(io::IO, sp::Species)
  out = string(gas(sp))
  open_bracket = false
  if !(charge(sp) isa Neutral) && string(gas(sp)) != "e"
    open_bracket = true
    out *= "[$(charge(sp))"
  end
  if !isnothing(electronic_state(sp))
    out *= !open_bracket ? "[" : ","
    open_bracket = true
    out *= string(electronic_state(sp))
  end
  if isnothing(vibrational_state(sp))
    out *= open_bracket ? "]" : ""
    return Base.print(io, out)
  else
    out *= ",vib=" * format_quantum_label(vibrational_state(sp))
  end
  if isnothing(rotational_state(sp))
    out *= open_bracket ? "]" : ""
    return Base.print(io, out)
  else
    out *= ",rot=" * format_quantum_label(rotational_state(sp)) * "]"
  end
  Base.print(io, out)
end

# Sorting

"""
    _quantum_label_isless(a, b)

Order two `QuantumLabel` values. Uses `isless` directly where Base defines it
(`Int`, `Tuple`, `Symbol`); falls back to comparing `string(a)`/`string(b)` for
pairs without a natural order (mixed types, `NamedTuple`) so `Species` sorting
stays total and never throws.
"""
function _quantum_label_isless(a, b)
  applicable(isless, a, b) && return isless(a, b)
  return isless(string(a), string(b))
end

Base.isless(::Negative, ::Charge) = false

Base.isless(::Positive, ::Negative) = true
Base.isless(::Neutral, ::Negative) = true

Base.isless(::Positive, ::Neutral) = false
Base.isless(::Neutral, ::Positive) = true

Base.isless(a::Negative, b::Negative) = isless(a.value, b.value)
Base.isless(a::Positive, b::Positive) = isless(a.value, b.value)

function Base.isless(a::Species, b::Species)
  if string(gas(a)) == "e"
    return false
  elseif string(gas(b)) == "e"
    return true
  end

  if charge(a) != charge(b)
    return isless(charge(a), charge(b))
  end

  if gas(a) != gas(b)
    return isless(gas(a), gas(b))
  end

  if electronic_state(a) != electronic_state(b)
    if electronic_state(a) == "X"
      return true
    end
    if electronic_state(b) == "X"
      return false
    end
    return isless(electronic_state(a), electronic_state(b))
  end

  if vibrational_state(a) != vibrational_state(b)
    if isnothing(vibrational_state(a))
      return false
    end
    if isnothing(vibrational_state(b))
      return true
    end
    return _quantum_label_isless(vibrational_state(a), vibrational_state(b))
  end

  if rotational_state(a) != rotational_state(b)
    if isnothing(rotational_state(a))
      return false
    end
    if isnothing(rotational_state(b))
      return true
    end
    return _quantum_label_isless(rotational_state(a), rotational_state(b))
  end
end
