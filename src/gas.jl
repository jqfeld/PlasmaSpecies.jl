abstract type Gas end

include("nuclide.jl")
include("nuclide_data.jl")
include("formula.jl")
include("gases/electron.jl")

"""
    REGISTERED_GASES

Labels that are not chemical formulas and so cannot be resolved by
[`parse_formula`](@ref). Only the electron qualifies.
"""
const REGISTERED_GASES = Dict{String,Any}(
  "e" => Electron,
)

"""
    Gas(label) -> Gas

Resolve a gas label, in three steps:

1. a special, non-formula label from [`REGISTERED_GASES`](@ref) — currently
   just `"e"`;
2. otherwise a chemical formula, via [`parse_formula`](@ref), giving a
   [`Molecule`](@ref) that knows its own mass;
3. otherwise [`StringGas`](@ref), so an unrecognised label still parses,
   prints and round-trips — it just cannot be weighed.
"""
function Gas(x)
  haskey(REGISTERED_GASES, x) && return REGISTERED_GASES[x]()
  parsed = parse_formula(x)
  isnothing(parsed) || return parsed
  return StringGas(x)
end

"""
    StringGas(name) <: Gas

Fallback for a label that is not a resolvable chemical formula. Behaves like any
other gas except that it has no composition, and therefore no mass.
"""
struct StringGas <: Gas
  name::String
end

Base.show(io::IO, g::StringGas) = print(io, g.name)

Base.isless(a::Gas, b::Gas) = isless(string(a), string(b))

composition(g::Gas) = error(
  """$(typeof(g)) has no composition. Only gases resolved from a chemical
  formula (see `parse_formula`) carry one.""")

mass(g::Gas) = error(
  """No mass available for $(g): its label does not resolve to a chemical formula.
  Check the spelling (symbols are case-sensitive, so `Co` is cobalt and `CO` is
  carbon monoxide), or add the missing nuclide to `ISOTOPE_MASSES` in
  src/nuclide_data.jl.""")
