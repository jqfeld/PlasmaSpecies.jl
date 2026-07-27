module PlasmaSpecies

import AbstractTrees

include("charge.jl")
export is_parent_species, get_parent_species
export Positive, Negative, Neutral, Charge, Species
export ispositive, isnegative, isneutral

include("gas.jl")
export Gas, StringGas, DiNitrogen, Nitrogen, DiOxygen, Oxygen

include("electronic_state.jl")
export ElectronicState, StringElectronicState, TermSymbol
export label

include("species.jl")
export gas, charge, electronic_state, vibrational_state, rotational_state, energy, degeneracy, mass
export QuantumLabel, level_matches, species_matches

include("species_tree.jl")
export SpeciesTree, leaves, apply_energy!, matching_leaves

include("reactions.jl")
export ReactionFormula, PlasmaReaction
export parse_reaction, @p_str
export apply_tree
export reaction_energy, thermal_source_term, isreactive


# CatalystExt definitions
to_catalyst(_...) = error("to_catalyst() not implemented. Is Catalyst.jl loaded?")
export to_catalyst

end 
