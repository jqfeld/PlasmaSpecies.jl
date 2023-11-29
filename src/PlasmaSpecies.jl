module PlasmaSpecies

import AbstractTrees

include("charge.jl")
export is_parent_species, get_parent_species
export Positive, Negative, Neutral, Charge, Species
export ispositive, isnegative, isneutral

include("gas.jl")
export Gas, StringGas, DiNitrogen, Nitrogen

include("electronic_state.jl")
export ElectronicState, StringElectronicState

include("species.jl")
export gas, charge, electronic_state, vibrational_state, rotational_state, mass

include("species_tree.jl")
export SpeciesTree, leaves

include("reactions.jl")
export PlasmaReaction
export parse_reaction, @p_str 
export apply_tree


# CatalystExt definitions
to_catalyst(_...) = error("to_catalyst() not implemented. Is Catalyst.jl loaded?")
export to_catalyst

end 
