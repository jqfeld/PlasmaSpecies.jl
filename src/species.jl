"""
    Species(
        gas <: Gas,
        charge=Neutral(),
        electronic_state=nothing,
        vibrational_state=nothing,
        rotational_state=nothing
    )

## Fields
#
- `gas <: Gas`: Label of the parent gas, e.g. a struct like Nitrogen <: Gas or the general StringGas(str:String).
- `charge=Neutral()`: Charge of the species, e.g. Neutral(), Positive(n), Negative(n)
- `electronic_state=nothing`: Optional, label for the electronic state.
- `vibrational_state=nothing`: Optional, label for the vibrational state. If defined, `electronic_state` cannot be `nothing`.
- `rotational_state=nothing`: Optional, label for the rotational state. If defined, `vibrational_state` cannot be `nothing`.
- `degeneracy=nothing`: Optional, statistical weight (degeneracy) of this state.

"""
Base.@kwdef struct Species{S,C<:Charge,E,V,J,En,D}
  gas::S
  charge::C = Neutral()
  electronic_state::E = nothing
  vibrational_state::V = nothing
  rotational_state::J = nothing
  energy::En = nothing
  degeneracy::D = nothing
end

"""
    QuantumLabel

The closed set of forms `vibrational_state`/`rotational_state` can take: a single
quantum number (`Int`), several quantum numbers (`Tuple`/`NamedTuple`, e.g. normal
modes of a polyatomic), a contiguous range of levels lumped together for reaction
purposes (`UnitRange{Int}`, see [`level_matches`](@ref)/[`matching_leaves`](@ref)),
or an opaque hand-picked label (`Symbol`) with no containment semantics.
"""
const QuantumLabel = Union{Int,NTuple{N,Int} where N,NamedTuple,UnitRange{Int},Symbol}

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


isnegative(x) = false
isnegative(::Species{S,Negative,E,V,J}) where {S,E,V,J} = true
ispositive(x) = false
ispositive(::Species{S,Positive,E,V,J}) where {S,E,V,J} = true
isneutral(x) = false
isneutral(::Species{S,Neutral,E,V,J}) where {S,E,V,J} = true

Base.:(==)(x::S1, y::S2) where {S1<:Species, S2<:Species} =
  gas(x) == gas(y) &&
  charge(x) == charge(y) &&
  electronic_state(x) == electronic_state(y) &&
  vibrational_state(x) == vibrational_state(y) &&
  rotational_state(x) == rotational_state(y)


"""
    Species(str::String)

Convenience constructor for the `Species` struct. 
It parses a string of the format defined by the [LoKI-B](https://github.com/IST-Lisbon/LoKI) 
Boltzmann solver and automatically fills in the fields of `Species`.
"""
function Species(str::String)
  str = join(map(x -> isspace(str[x]) ? "" : str[x], 1:length(str)))
  regex = r"(^(?:\w|\d)*)(?:\[([+\-]+)?,?([^,\]]+)?(?:,vib=(\([^\)]*\)|[^,\]]*)(?:,rot=(\([^\)]*\)|[^,\]]*))?)?\])?"
  m = collect(match(regex, str))
  Species(
    gas = Gas(m[1]),
    charge = m[1] == "e" ? Negative(1) : Charge(m[2]),
    electronic_state = ElectronicState(m[3]),
    vibrational_state = parse_quantum_label(m[4]),
    rotational_state = parse_quantum_label(m[5]))
    # m[4:end]...)
end
Species(gas::G) where {G<:Gas} = Species(gas, Neutral(), nothing, nothing, nothing, nothing, nothing)


gas(sp::Species) = sp.gas
charge(sp::Species) = sp.charge
electronic_state(sp::Species) = sp.electronic_state
vibrational_state(sp::Species) = sp.vibrational_state
rotational_state(sp::Species) = sp.rotational_state
energy(sp::Species) = sp.energy
degeneracy(sp::Species) = sp.degeneracy
mass(sp::Species) = gas(sp) isa Electron ? mass(gas(sp)) : mass(gas(sp)) - to_value(charge(sp)) * mass(Electron())

function get_parent_species(sp::Species)
  if !isnothing(rotational_state(sp))
    return Species(
      gas(sp),
      charge(sp),
      electronic_state(sp),
      vibrational_state(sp),
      nothing,
      nothing,
      nothing
    )
  elseif !isnothing(vibrational_state(sp))
    return Species(
      gas(sp),
      charge(sp),
      electronic_state(sp),
      nothing,
      nothing,
      nothing,
      nothing
    )
  elseif !isnothing(electronic_state(sp))
    return Species(
      gas(sp),
      charge(sp),
      nothing,
      nothing,
      nothing,
      nothing,
      nothing
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

function is_parent_species(sp::Species, parent::Species)
  if parent != sp && gas(parent) == gas(sp)
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
