module PlasmaSpecies

import AbstractTrees


include("species.jl")
export is_parent_species, get_parent_species
export Positive, Negative, Neutral, Charge, Species
export ispositive, isnegative, isneutral
export gas, charge, electronic_state, vibrational_state, rotational_state

include("species_tree.jl")
export SpeciesTree, leaves

include("reactions.jl")
export PlasmaReaction
export parse_reaction, @p_str 
export apply_tree


# CatalystExt definitions
to_catalyst(_...) = error("This function needs Catalyst.jl")
export to_catalyst

end 
