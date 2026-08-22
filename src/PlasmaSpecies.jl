module PlasmaSpecies

import AbstractTrees

include("charge.jl")
export is_parent_species, get_parent_species
export Positive, Negative, Neutral, Charge, Species
export ispositive, isnegative, isneutral

include("gas.jl")
export Gas, StringGas, Electron, Molecule
export Nuclide, Element, Isotope, ELEMENTS, ISOTOPES
export parse_formula, composition, canonical_composition
export element, isotope, atomic_number, mass_number

include("electronic_state.jl")
export ElectronicState, StringElectronicState, TermSymbol
export label

include("species.jl")
export gas, charge, electronic_state, vibrational_state, rotational_state, energy, degeneracy, metadata, mass
export with_fields
export QuantumLabel, level_matches, species_matches

include("species_tree.jl")
export SpeciesTree, leaves, apply_energy!, matching_leaves

include("reactions.jl")
export ReactionFormula, PlasmaReaction
export parse_reaction, @p_str
export apply_tree
export reaction_energy, thermal_source_term, isreactive, scale_rate


"""
    to_catalyst(sp::Species)
    to_catalyst(t::SpeciesTree)
    to_catalyst(t::SpeciesTree, reactions...)

Export to [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/). Provided by an
extension: it only works once `Catalyst` is loaded in the same session,
otherwise every method throws.

- A `Species` becomes a Catalyst species variable named after its string, so
  anything carrying a charge or a state label is a `var"..."` identifier
  (`var"N2[+]"`, `var"N2[X,vib=1]"`).
- A `SpeciesTree` becomes the sum of its leaf variables.
- `reactions` are `(rate, formula)` tuples or [`PlasmaReaction`](@ref)s, given
  as varargs or a vector. Each is expanded over the tree by [`apply_tree`](@ref)
  — branching factors folded into the rate — and returned as a
  `Vector{Reaction}` ready for a `ReactionSystem`. Reactions that throw are
  skipped with a warning.
"""
to_catalyst(_...) = error("to_catalyst() not implemented. Is Catalyst.jl loaded?")
export to_catalyst

end 
